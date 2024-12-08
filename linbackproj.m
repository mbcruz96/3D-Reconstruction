function [features, valid_points] = getFeatures(img)
%{
    This function takes an image and finds the 200 strongest features
    using SIFT. Using those points, each points feature vector and location
    are extracted and returned.
%}

% detecting all SIFT features
points = detectSURFFeatures(rgb2gray(img));
% selecting >= 1000 features from the points
points = points.selectStrongest(1000);
% extracting the feature vectors and locations from strongest points
[features, valid_points] = extractFeatures(rgb2gray(img), points);
end

function [match1, match2] = getMatchingPoints(features1, points1, features2, points2)
%{
    This function takes the feature vectors and locations from the
    strongest SIFT features from two images and finds matching features 
    between the images. It returns two matrices corresponding to the 
    homogeneous coordinates of the matching points in a two vector. 
%}

% finding matching features between both images and storing their indexes
index_pairs = matchFeatures(features1, features2);

% extracting matching points
match1 = points1(index_pairs(:, 1), :);
match2 = points2(index_pairs(:, 2), :);

% making each point homogenous
match1 = [match1.Location, ones(size((match1.Location), 1), 1)];
match2 = [match2.Location, ones(size((match2.Location), 1), 1)];
end

function correspondence_matrix = constructCorrespondenceMatrix(x, x_prime)
%{
    This function constructs the point correspondence matrix. The function
    accepts two matrices of matching points x and x_prime and for each set
    of matching points constructs the 1x9 vector:
    [x'x,  x'y,  x',  y'x,  y'y,  y',  x,  y,  1]
    The vectors for each pair of matching points are stacked ontop of
    eachother, which constructs the correspondence matrix
%}

% matrix with all point correspondences
correspondence_matrix = single.empty;

% iterating through all matching points
for i = 1:size(x, 1)
    correspondence_matrix = [correspondence_matrix; x_prime(i, 1)*x(i, 1), x_prime(i, 1)*x(i,2), x_prime(i, 1), x_prime(i, 2)*x(i, 1), x_prime(i, 2)*x(i, 2), x_prime(i, 2), x(i, 1), x(i, 2), 1];
end

end

function fundamental_matrix = GetFundamentalMatrix(correspondence_matrix)
%{
    This function takes in a correspondence matrix representing the A
    matrix in the equation Ah = 0. To find the nullspace of A, SVD will be 
    performed on matrix A and the last column of the V matrix is the
    solution and will be returned. 
%}

% performing SVD on A matrix
[U, D, V] = svd(correspondence_matrix);
% extracting last column of V matrix and constructing 3x3 matrix
fundamental_matrix = [V(1:3, 9)'; V(4:6, 9)'; V(7:9, 9)'];
[U, D, V] = svd(fundamental_matrix);
% enforcing singularity constraint
D(3, 3) = 0;
% constructing singular fundamental matrix
fundamental_matrix = U * D * V.';

end

function fundamental_matrix = EstimateFundamentalMatrix(match1, match2)
%{
    This function estimates the best homography for a set of matching
    points by doing a modified version of RANSAC and iteratively finding inliers between the two matching point
    sets and extracting the coordinates for the inliers. Inliers are
    extracted until there are only 4 inliers left, or the number of inliers
    doesn't change between iterations. Once complete, the set of points is
    used to calculate the homography at infinity between the images which is returned. 
%}
m1 = match1;
m2 = match2;
num_points = size(m1, 1); % current total number of inliers
prev_points = -1; % previous number of inliers

% select inliers until the amount of inliers is 4 or the number of inliers
% doesn't change
while (num_points > 8) && (prev_points ~= num_points)
    % finding inliers
    [~, inlier_idx] = estgeotform2d(m1(:, 1:2), m2(:, 1:2), "projective");
    % extracting the indexes of inliers
    indexes = find(inlier_idx == 1);
    % extracting the coordinates of the inliers
    m1 = m1(indexes, : );
    m2 = m2(indexes, : );
    % adjusting stop critera variables
    prev_points = num_points;
    num_points = size(m1, 1);
end

% constructing the fundamental matrix from inliers
correspondence = constructCorrespondenceMatrix(m1, m2);
fundamental_matrix = GetFundamentalMatrix(correspondence);
end

function [P1, P2] = GetCameraMatrix(F, intrinsics, match1, match2) 
%{
    This function calculates the camera calibration matrices P1 and P2 for two
    cameras by constructing the matrices P1 = K[I|0] and P2 = K[R|t] where
    K is the camera's intrinsic parameters, R is the rotation matrix of
    camera 2, t is the translation matrix of camera 2, and I is the
    identity matrix. The function accepts the fundamental matrix between
    the two scenes F, the camera intrinsic parameters K, and a set of
    matching points from camera 1 match1 and camera 2 match2. 
%}

% Getting orientation and location of camera 2 in reference to camera 1
[relativeOrientation, relativeLocation] = relativeCameraPose(F, intrinsics, match1(:, 1:2), match2(:, 1:2));
% Getting rotation matrix and translation vector
[rotationMatrix, translationMatrix] = cameraPoseToExtrinsics(relativeOrientation, relativeLocation);

% Constructing [I|0] matrix
I = [eye(3), zeros(3,1)];
% Constructing [R|t] matrix
R_t = [rotationMatrix, translationMatrix.'];

% Calculating camera calibration matrix
P1 = intrinsics.K*I;
P2 = intrinsics.K*R_t;

end

function points = get3DPoints(x, x_prime, p, p_prime)
%{
    This function will compute 3D coordinates from two matching points x
    and x_prime and two camera projection matricies p and p_prime by 
    constructing matrices x x PX = 0 and x' x P'X = 0 and stacking them 
    together. Then the singular value decomposition of the matrix will be 
    computed to get the null space representing the 3D coordinate.
%}

% initializing correspondence matrix
points = single.empty;

% looping through point correspondences
for i = 1:size(x, 1)
    % calculating 3x9 submatrix for each point and concatenating points
    point_matrix = [x(i, 1)*p(3, :) - p(1, :)
                    x(i, 2)*p(3, :) - p(2, :)
                    %x(i, 1)*p(2, :) - p(1, :)
                    x_prime(i, 1)*p_prime(3, :) - p_prime(1, :)
                    x_prime(i, 2)*p_prime(3, :) - p_prime(2, :)];
                    %x_prime(i, 1)*p_prime(2, :) - p_prime(1, :)];
    [U, D, V] = svd(point_matrix);
    point = V(1:3, 4).'  / V(4,4);
    points = [points; point];
end

end

% running calibration script for specified image folder
image_dir = 'Playroom';
run(strcat('Assignment4/', image_dir, '/calib.m'));


% Getting feature points for each image
[features1, points1] = getFeatures(im1);
[features2, points2] = getFeatures(im2);

% Getting matching points of the two images
[match1, match2] = getMatchingPoints(features1, points1, features2, points2);
figure; 
showMatchedFeatures(im1, im2, match1(:, 1:2), match2(:, 1:2));

% Estimating the fundamental matrix
fundamental_matrix = EstimateFundamentalMatrix(match1, match2);

% Constructing camera intrinsic parameter matrix
focal_length = [cam2(1,1), cam2(2,2)];
principal_point = [cam2(1, 3), cam2(2, 3)];
image_size = [size(im1, 1), size(im1, 2)];
intrinsic_params = cameraIntrinsics(focal_length, principal_point, image_size);

% Constructing camera projection matrices
[P1, P2] = GetCameraMatrix(fundamental_matrix, intrinsic_params, match1, match2);

% Getting 3D coordinates for all matching points
point_cloud = get3DPoints(match1(:, 1:2), match2(:, 1:2), P1, P2);

% Computing the point cloud model 
pcshow(point_cloud, [0 1 0]);


% verifying results with triangulate 
world_points = triangulate(match1(:, 1:2), match2(:, 1:2), P1, P2);
%pcshow(world_points, [0 1 0]);
