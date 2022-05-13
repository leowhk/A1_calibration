/**
 * Copyright (C) 2015 by Liangliang Nan (liangliang.nan@gmail.com)
 * https://3d.bk.tudelft.nl/liangliang/
 *
 * This file is part of Easy3D. If it is useful in your research/work,
 * I would be grateful if you show your appreciation by citing it:
 * ------------------------------------------------------------------
 *      Liangliang Nan.
 *      Easy3D: a lightweight, easy-to-use, and efficient C++
 *      library for processing and rendering 3D data. 2018.
 * ------------------------------------------------------------------
 * Easy3D is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License Version 3
 * as published by the Free Software Foundation.
 *
 * Easy3D is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */


#include "calibration.h"
#include "matrix_algo.h"

using namespace easy3d;

/*
 * TODO: Finish this function for calibrating a camera from the corresponding 3D-2D point pairs.
 *       You may define a few functions for some sub-tasks.
 * @return True on success, otherwise false. On success, the camera parameters are returned by
 */

Matrix construct_P(const std::vector<Vector3D>& pnt_3d,
                   const std::vector<Vector2D>& pnt_2d){
    Matrix P_mat( 2 * pnt_3d.size(), 12,0.0);
    std::vector<double> P_xi;
    std::vector<double> P_yi;
    for (int i=0; i < pnt_3d.size() ; i++){
        double X_3d = pnt_3d[i].x();
        double Y_3d = pnt_3d[i].y();
        double Z_3d = pnt_3d[i].z();
        double x_img = pnt_2d[i].x();
        double y_img = pnt_2d[i].y();

        P_xi = {X_3d, Y_3d, Z_3d, 1, 0, 0, 0, 0, -x_img * X_3d, -x_img * Y_3d, -x_img * Z_3d, -x_img};
        P_yi = {0, 0, 0, 0, X_3d, Y_3d, Z_3d, 1, -y_img * X_3d, -y_img * Y_3d, -y_img * Z_3d, -y_img};

        std::cout << P_xi << std::endl;
        std::cout << P_yi << std::endl;

        P_mat.set_row(2*i, P_xi);
        P_mat.set_row(2*i+1, P_yi);
    }
    std::cout << P_mat << std::endl;
    return P_mat;
}


bool Calibration::calibration(
        const std::vector<Vector3D>& points_3d, /// input: An array of 3D points.
        const std::vector<Vector2D>& points_2d, /// input: An array of 2D image points.
        double& fx, double& fy,    /// output: the focal length (in our slides, we use 'alpha' and 'beta'),
        double& cx, double& cy,    /// output: the principal point (in our slides, we use 'u0' and 'v0'),
        double& skew,              /// output: the skew factor ('-alpha * cot_theta')
        Matrix33& R,               /// output: the 3x3 rotation matrix encoding camera orientation.
        Vector3D& t)               /// output：a 3D vector encoding camera translation.
{

    // TODO: check if input is valid (e.g., number of correspondences >= 6, sizes of 2D/3D points must match)
    if (points_3d.size() >= 6 || points_2d.size() >= 6) {

        // TODO: construct the P matrix (so P * m = 0).

        Matrix P_mat = construct_P(points_3d, points_2d);

        // TODO: solve for M (the whole projection matrix, i.e., M = K * [R, t]) using SVD decomposition.
        //   Optional: you can check if your M is correct by applying M on the 3D points. If correct, the projected point
        //             should be very close to your input images points.

        Matrix U_mat(P_mat.rows(), P_mat.rows(), 0.0);
        Matrix S_mat(P_mat.rows(), 12, 0.0);
        Matrix V_mat(12, 12, 0.0);

        std::cout << "SVD decomposition ...." << std::endl;
        svd_decompose(P_mat, U_mat, S_mat, V_mat);
        std::cout << U_mat << std::endl;
        std::cout << S_mat << std::endl;
        std::cout << V_mat << std::endl;

        std::cout << "last column of V .... is M" << std::endl;
        Vector M_vec = V_mat.get_column(V_mat.cols() - 1);
        std::cout << M_vec << std::endl;

        // TODO: extract intrinsic parameters from M.
        std::cout << "M_mat" << std::endl;
        Matrix M_mat(3, M_vec.size() / 3, M_vec.data());
        std::cout << M_mat << std::endl;

        Matrix A_mat(3, 3, 0.0);
        for (int i = 0; i < M_mat.cols() - 1; i++) {
            Vector A_vector = M_mat.get_column(i);
            A_mat.set_column(i, A_vector);
        }

        std::cout << "A_mat" << std::endl;
        std::cout << A_mat << std::endl;

        Vector A1_vec = A_mat.get_row(0);
        Vector A2_vec = A_mat.get_row(1);
        Vector A3_vec = A_mat.get_row(2);

        std::cout << "A1_vec" << std::endl;
        std::cout << A1_vec << std::endl;
        std::cout << A2_vec << std::endl;
        std::cout << A3_vec << std::endl;


        double po = 1 / norm(A3_vec);

        cx = pow(po, 2) * dot(A1_vec, A2_vec);
        cy = pow(po, 2) * dot(A2_vec, A3_vec);

        std::cout << "po: " << po << std::endl;
        std::cout << "cx: " << cx << std::endl;
        std::cout << "xy: " << cy << std::endl;

        double nom = dot(cross(A1_vec, A3_vec), cross(A2_vec, A3_vec));
        double denom = norm(cross(A1_vec, A3_vec)) * norm(cross(A2_vec, A3_vec));

        double theta = acos(-1 * nom / denom);

        fx = pow(po, 2) * norm(cross(A1_vec, A2_vec)) * sin(theta);
        fy = pow(po, 2) * norm(cross(A2_vec, A3_vec)) * sin(theta);
        skew = -1 * fx / tan(theta);

        std::cout << "theta: " << theta << std::endl;
        std::cout << "skew: " << skew << std::endl;

        // TODO: extract extrinsic parameters from M.

        Vector R3_vec = po * A3_vec;
        Vector R1_vec = cross(A2_vec, A3_vec) / norm(cross(A2_vec, A3_vec));
        Vector R2_vec = cross(R3_vec, R1_vec);

        R.set_row(0, R1_vec);
        R.set_row(1, R2_vec);
        R.set_row(2, R3_vec);

        double cotan = cos(theta) / sin(theta);

        Matrix33 K = {fx, -fx * cotan, cx,
                      0, fx / sin(theta), cx,
                      0, 0, 1};

        Vector b_mat = M_mat.get_column(M_mat.cols() - 1);

        t = po * inverse(K) * b_mat;

        std::cout << R << std::endl;
        std::cout << t << std::endl;

        std::cout
                << "\n\tTODO: After you implement this function, please return 'true' - this will trigger the viewer to\n"
                   "\t\tupdate the rendering using your recovered camera parameters. This can help you to visually check\n"
                   "\t\tif your calibration is successful or not.\n\n" << std::flush;
        return true;
    }

    else{
        return false;
    }
}


