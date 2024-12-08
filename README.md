# Stereo 3D-Reconstruction
- Stereo 3D Reconstruction by linear triangulation implementation
- For each image file, pairs of matching points are found between the two image views of the scene
- Using a modified version of RANSAC, inliers of matching points are detected
- The fundamental matrix between the scenes is computed using the inlier matching image points
- The camera project matrix is constructed for each image using the fundamental matrix, camera intrinsic parameters defined in the calib.m file of each image folder, and the set of matching points
- The 3D coordinate of each pair of matching image points is computed by a linear triangulation algorithm
- The collection of 3D points are stacked into a point cloud model and the 3D reconstruction is rendered
- Another 3D reconstruction of the scene is created with the built-in MatLab Triangulation function for comparison
- All results can be viewed in the experimental report and the results directory
