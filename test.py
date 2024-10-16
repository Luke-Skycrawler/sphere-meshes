from sqem import Sqem
import numpy as np

points = np.array(
    [[1.0, 0.0, 0.0], 
    [-1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0], 
    [0.0, -1.0, 0.0],
    [0.0, 0.0, 1.0], 
    [0.0, 0.0, -1.0]])


normals = points
# cubes

sqems = [Sqem(points[i], normals[i]) for i in range(points.shape[0])]
sum = sqems[0]
for i in range(1, points.shape[0]):
    sum = sum + sqems[i]

center = np.ones(3)
radius = float(-1.0)
pa = np.ones(3) * -10.0
pb = np.ones(3) * 10.0


center, radius = sum.minimize(pa, pb)
print(center, radius)

points2 = np.array(
    [[2.5, 2., 0.],
     [1.5, 2., 0.],
     [2., 1.5, 0.],
     [2., 2., 0.5],
     [2., 2., -0.5]])

sqems2 = [Sqem(points2[i], normals[i]) for i in range(points2.shape[0])]
sum2 = sqems2[0]
for i in range(1, points2.shape[0]):
    sum2 = sum2 + sqems2[i]

center2, radius2 = sum2.minimize(pa, pb)
print(center2, radius2)