import numpy as np
import scipy.spatial as ss

# Face 0
points = np.array([[0, 0, 0],
                    [0.05, 0, 0],   # here
                    [0.1, 0, 0],
                    [0, 0.05, 0],
                    [0.05, 0.05, 0],    # here
                    [0.1, 0.05, 0],
                    [0, 0.1, 0],
                    [0.05, 0.1, 0],
                    [0.1, 0.1, 0],
                    [0, 0, 0.01],
                    [0.05, 0, 0.01], # here
                    [0.1, 0, 0.01],
                    [0, 0.05, 0.01],
                    [0.05, 0.05, 0.01], #here
                    [0.1, 0.05, 0.01],
                    [0, 0.1, 0.01],
                    [0.05, 0.1, 0.01],
                    [0.1, 0.1, 0.01]])
face_points_ids = np.array([[1, 4, 13, 10],
                            [3, 12, 13, 4],
                            [4, 13, 14, 5],
                            [4, 7, 16, 13],
                            [6, 15, 16, 7],
                            [7, 16, 17, 8],
                            [0, 9, 12, 3],
                            [3, 12, 15, 6],
                            [2, 5, 14, 11],
                            [5, 8, 17, 14],
                            [0, 1, 10, 9],
                            [1, 2, 11, 10],
                            [0, 3, 4, 1],
                            [3, 6, 7, 4],
                            [1, 4, 5, 2],
                            [4, 7, 8, 5],
                            [9, 10, 13, 12],
                            [12, 13, 16, 15],
                            [10, 11, 14, 13],
                            [13, 14, 17, 16]])

face_zero_points = [1, 4, 13, 10]

centroid_face = np.array([0., 0., 0.])
for pi in face_zero_points:
    centroid_face += points[pi]
centroid_face /= len(face_zero_points)

print(centroid_face)

# Face area checked (and corrected) manually
cell_1_face_ids = np.array([0, 1, 6, 10, 12, 16])

all_cell_points = []
for f_id in cell_1_face_ids:
    all_cell_points += list(face_points_ids[f_id])

new_list = []
for i in all_cell_points:
    if i not in new_list:
        new_list.append(i)

print(new_list)
centroid_cell = np.array([0., 0., 0.])
new_list_points = np.array([0., 0., 0.])
for pi in new_list:
    centroid_cell += points[pi]
    new_list_points = np.vstack((new_list_points, points[pi]))
centroid_cell /= len(new_list)
print(centroid_cell)

hull = ss.ConvexHull(new_list_points[1:, :])
print('volume inside points is: ',hull.volume)