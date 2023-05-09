import sys
import math
import os
import numpy as np
import scipy.spatial as spatial


def make_dir_for_file(file):
    file_dir = os.path.dirname(file)
    if file_dir != '':
        if not os.path.exists(file_dir):
            try:
                os.makedirs(os.path.dirname(file))
            except OSError as exc: # Guard against race condition
                raise

def write_xyz(file_path, points: np.ndarray, normals=None, colors=None):
    """
    Write point cloud file.
    :param file_path:
    :param points:
    :param normals:
    :param colors:
    :return: None
    """

    make_dir_for_file(file_path)

    if points.shape == (3,):
        points = np.expand_dims(points, axis=0)

    if points.shape[0] == 3 and points.shape[1] != 3:
        points = points.transpose([1, 0])

    if colors is not None and colors.shape[0] == 3 and colors.shape[1] != 3:
        colors = colors.transpose([1, 0])

    if normals is not None and normals.shape[0] == 3 and normals.shape[1] != 3:
        normals = normals.transpose([1, 0])

    with open(file_path, 'w') as fp:

        # convert 2d points to 3d
        if points.shape[1] == 2:
            vertices_2p5d = np.zeros((points.shape[0], 3))
            vertices_2p5d[:, :2] = points
            vertices_2p5d[:, 2] = 0.0
            points = vertices_2p5d

        # write points
        # meshlab doesn't like colors, only using normals. try cloud compare instead.
        for vi, v in enumerate(points):
            line_vertex = str(v[0]) + " " + str(v[1]) + " " + str(v[2]) + " "
            if normals is not None:
                line_vertex += str(normals[vi][0]) + " " + str(normals[vi][1]) + " " + str(normals[vi][2]) + " "
            if colors is not None:
                line_vertex += str(colors[vi][0]) + " " + str(colors[vi][1]) + " " + str(colors[vi][2]) + " "
            fp.write(line_vertex + "\n")


def load_xyz(file_path):
    data = np.loadtxt(file_path).astype('float32')
    nan_lines = np.isnan(data).any(axis=1)
    num_nan_lines = np.sum(nan_lines)
    if num_nan_lines > 0:
        data = data[~nan_lines]  # filter rows with nan values
        print('Ignored {} points containing NaN coordinates in point cloud {}'.format(num_nan_lines, file_path))
    return data



if __name__ == '__main__':
    
    input_path = "./in"
    out_path = "./out"
    all_xyzs = os.listdir(input_path)
    expected_num = 60000
    for xyz in all_xyzs:
        in_pc_path = os.path.join(input_path, xyz)
        out_pc_path = os.path.join(out_path, xyz)
        pc = load_xyz(in_pc_path)
        point_num = pc.shape[0]
        if point_num <= expected_num:
            print(xyz)
            write_xyz(out_pc_path, points=pc[:,:3], normals=pc[:,3:])
            continue
        
        
        sub_index = np.random.choice(point_num, expected_num, replace=False) #randomly subsample a point cloud with expected_num points
        pc_new = pc[sub_index, :]
        
        
        print(xyz)
        write_xyz(out_pc_path, points=pc_new[:,:3], normals=pc_new[:,3:])

