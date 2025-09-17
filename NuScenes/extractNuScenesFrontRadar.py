# this is only front radar
import os
import numpy as np
from scipy.io import savemat
from nuscenes.nuscenes import NuScenes
from nuscenes.utils.data_classes import RadarPointCloud
from nuscenes.map_expansion.map_api import NuScenesMap
from PIL import Image

# === Config ===
DATAROOT = 'D://nuscenesData/v1.0-mini'
VERSION = 'v1.0-mini'
OUT_FOLDER = 'nuscenes_scene_mat_output'
os.makedirs(OUT_FOLDER, exist_ok=True)

# === Load NuScenes ===
nusc = NuScenes(version=VERSION, dataroot=DATAROOT, verbose=True)

# === Loop through each scene ===
for scene in nusc.scene:
    scene_name = scene['name']
    log = nusc.get('log', scene['log_token'])

    # nuScenes older minis use 'location'; newer use 'map_name'
    location_to_map = {
        'boston-seaport': 'boston-seaport',
        'singapore-onenorth': 'singapore-onenorth',
        'singapore-hollandvillage': 'singapore-hollandvillage',
        'singapore-queenstown': 'singapore-queenstown'
    }
    map_name = log.get('map_name') or location_to_map.get(log.get('location'))
    nusc_map = NuScenesMap(dataroot=DATAROOT, map_name=map_name)

    print(f"\nProcessing scene: {scene_name}")
    sample_token = scene['first_sample_token']
    scene_samples = []

    while sample_token:
        sample = nusc.get('sample', sample_token)

        try:
            # === Radar & Camera tokens ===
            radar_token = sample['data']['RADAR_FRONT']
            cam_token   = sample['data']['CAM_FRONT']

            # === Radar Front Point Cloud (single sweep for 1:1 parity with MATLAB) ===
            radar_sd = nusc.get('sample_data', radar_token)
            radar_pc = RadarPointCloud.from_file(nusc.get_sample_data_path(radar_token))

            # radar_pc.points: 18 x N (channels x points)
            p = radar_pc.points.T.astype(np.float32)  # N x 18

            # Columns (0-based): 0:x,1:y,2:z, 8:vx_comp, 9:vy_comp
            x  = p[:, 0]
            y  = p[:, 1]
            z  = p[:, 2]
            vx = p[:, 8]  # ego-motion–compensated vx
            vy = p[:, 9]  # ego-motion–compensated vy
            vz = np.zeros_like(vx, dtype=np.float32)  # radar is 2D velocity; set vz = 0

            # Stack to N x 6: [x, y, z, vx, vy, vz]
            radar_points = np.stack([x, y, z, vx, vy, vz], axis=1).astype(np.float32)

            # === Radar Calibration (sensor → ego frame) ===
            radar_calib = nusc.get('calibrated_sensor', radar_sd['calibrated_sensor_token'])
            radar_translation = np.array(radar_calib['translation'], dtype=np.float32)  # radar → ego
            radar_rotation    = np.array(radar_calib['rotation'], dtype=np.float32)     # quaternion [w, x, y, z]

            # === Front Camera ===
            image_path = nusc.get_sample_data_path(cam_token)
            image_np = np.array(Image.open(image_path).convert('RGB'))

            # === Ego Pose ===
            ego_pose_token = radar_sd['ego_pose_token']
            ego_pose = nusc.get('ego_pose', ego_pose_token)
            ego_translation = np.array(ego_pose['translation'], dtype=np.float32)
            ego_rotation    = np.array(ego_pose['rotation'], dtype=np.float32)

            # === Map Lanes near ego ===
            ex, ey = float(ego_translation[0]), float(ego_translation[1])
            lane_dict = nusc_map.get_records_in_radius(ex, ey, 50, ['lane'])
            lanes = []
            for lane_token in lane_dict['lane']:
                lane = nusc_map.get('lane', lane_token)
                poly = nusc_map.extract_polygon(lane['polygon_token'])
                coords = np.array(poly.exterior.coords, dtype=np.float32)
                lanes.append(coords)
            lanes = np.array(lanes, dtype=object)

            # === Ground Truth Boxes ===
            gt_boxes = []
            for ann_token in sample['anns']:
                ann = nusc.get('sample_annotation', ann_token)
                gt_boxes.append({
                    'center':   np.array(ann['translation'], dtype=np.float32),
                    'size':     np.array(ann['size'], dtype=np.float32),
                    'rotation': np.array(ann['rotation'], dtype=np.float32),
                    'category': ann['category_name'],
                    'instance_token': ann['instance_token']
                })
            gt_boxes = np.array(gt_boxes, dtype=object)

            # === Store Sample ===
            sample_data = {
                # Save as N x 6: x, y, z, vx, vy, vz (vz=0)
                'radar_points': radar_points,

                # Poses / transforms
                'ego_translation': ego_translation,
                'ego_rotation':    ego_rotation,
                'radar_translation': radar_translation,
                'radar_rotation':    radar_rotation,

                # Map + image + bookkeeping
                'lanes': lanes,
                'camera_image_data': image_np,
                'sample_token': sample['token'],
                'gt_boxes': gt_boxes,
                'timestamp': radar_sd['timestamp'],

                'scene_description': scene.get('description', '')
                
            }
            #print(scene.get('description', ''))

            scene_samples.append(sample_data)

        except Exception as e:
            print(f" Error in sample {sample['token']}: {e}")

        sample_token = sample['next']

    # === Save Scene ===
    out_path = os.path.join(OUT_FOLDER, f"{scene_name}.mat")
    savemat(out_path, {'samples': np.array(scene_samples, dtype=object)}, do_compression=True)
    print(f"Scene saved: {out_path}")