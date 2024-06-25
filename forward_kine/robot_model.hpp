#ifndef ROBOT_MODEL_HPP_
#define ROBOT_MODEL_HPP_
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Geometry>
#include <list>
#include <cstdint>
#include <memory>

/*

        ↑進行方向

     _ _ _ _ _ _ _
    |             |
    |             |       　　右手
    |             |
    |_____________|

    ↑ +x
    ◎ → +y
    画面手前側が+z
    座標原点はロボットの左右腰ヨー軸の中心点

    単位 = mm
    座標原点のリンク番号 = 0
*/

constexpr double L01_Y_WAIST_YAW = 50.0;                  // 原点から腰ヨー軸までのy 腰ヨー軸間の距離 = * 2  最早いらないかも
constexpr double L12_Z_WAIST_YAW_TO_PITCH = 50.0;         // 腰ヨー軸から腰付け根ユニットまでのz
constexpr double L34_Z_WAIST_PITCH_TO_UPPER_KNEE = 132.5; // 腰付け根ユニットから膝上までのz
constexpr double L45_Z_UPPER_KNEE_TO_LOWER_KNEE = 55.0;   // 膝上から膝下までのz
constexpr double L56_Z_LOWER_KNEE_TO_FOOT_PITCH = 132.5;  // 膝下から足首までのz
constexpr double L67_Z_FOOT_PITCH_TO_FOOT_BOARD = 54.0;   // 足首から足裏までのz

/**
 * @brief ZYXオイラー角から回転行列を生成する
 * @param rotate ZYXオイラー角[rad]。渡すのはx,y,zの順番で良い。内部でZYXの順に掛ける
 * @return Eigen::Matrix3d 回転行列
 */
Eigen::Matrix4d get3DRotationMatrix(const Eigen::Vector3d &rotate)
{
    Eigen::Matrix4d rotate_matrix = Eigen::Matrix4d::Identity();
    Eigen::Matrix3d tmp;
    tmp = Eigen::AngleAxisd(rotate[2], Eigen::Vector3d::UnitZ()) *
          Eigen::AngleAxisd(rotate[1], Eigen::Vector3d::UnitY()) *
          Eigen::AngleAxisd(rotate[0], Eigen::Vector3d::UnitX());
    rotate_matrix.block<3, 3>(0, 0) = tmp;
    return rotate_matrix;
}

constexpr long double operator"" _deg(long double deg)
{
    return deg * (M_PI / 180.0);
}

struct StraightChainRobotModel
{
    static const uint32_t NO_PARALLEL_ID = -1;
    uint32_t motor_id_;
    uint32_t id_parallel_with_; // 平行リンクの対となるID
    bool is_dependent_;         // 平行リンクの場合、それが従属側か否か
    // DH記法で表現.
    Eigen::Vector3d x_link_length_; // x,y,z [mm]
    Eigen::Matrix4d x_link_rotate_; // ZYXオイラー角 [rad]
    Eigen::Vector3d z_link_length_;
    Eigen::Matrix4d z_link_rotate_;
    double joint_angle_; // [rad]
    bool does_reverse_;  // 1 or -1
    bool is_imagenary_;  // 動かさない仮想のリンク
    StraightChainRobotModel(const uint32_t motor_id, const bool does_reverse) : motor_id_(motor_id), does_reverse_(does_reverse_)
    {
    }
    StraightChainRobotModel(
        const uint32_t motor_id,
        const double x_link_length,
        const Eigen::Vector3d x_link_rotate,
        const double z_link_length,
        const Eigen::Vector3d z_link_rotate) : motor_id_(motor_id), id_parallel_with_(NO_PARALLEL_ID), is_dependent_(false),
                                               x_link_length_(x_link_length, 0, 0),
                                               x_link_rotate_(get3DRotationMatrix(x_link_rotate)),
                                               z_link_length_(0, 0, z_link_length),
                                               z_link_rotate_(get3DRotationMatrix(z_link_rotate)),
                                               is_imagenary_(false), does_reverse_(false)
    {
    }
    void printLinkState()
    {
        std::cout << "------- motor id <" << motor_id_ << "> --------" << std::endl;
        std::cout << "------- reverse <" << std::boolalpha << does_reverse_ << "> --------" << std::endl;
        std::cout << "*** x_link_length ***\n"
                  << x_link_length_ << std::endl;
        std::cout << "*** x_link_rotate ***\n"
                  << x_link_rotate_ << std::endl;
        std::cout << "*** z_link_length ***\n"
                  << z_link_length_ << std::endl;
        std::cout << "*** z_link_rotate ***\n"
                  << z_link_rotate_ << std::endl;
        std::cout << "-------------------------------------------------" << std::endl;
    }

    void setYoffset(const double y_offset)
    {
        x_link_length_(1, 0) = y_offset;
        return;
    }

    void setReverse(const bool does_reverse)
    {
        does_reverse_ = does_reverse;
        return;
    }

    void setImagenary(const bool is_imagenary)
    {
        is_imagenary_ = is_imagenary;
        return;
    }

    bool isParallel()
    {
        return id_parallel_with_ != NO_PARALLEL_ID;
    }

    bool isDependent()
    {
        return is_dependent_;
    }

    /**
     * @brief Set the Parallel With object
     *
     * @param id_parallel_with 平行リンクの相手のID
     * @param is_dependent 自分が従属側か否か
     */
    void setParallelWith(const uint32_t id_parallel_with, const bool is_dependent)
    {
        id_parallel_with_ = id_parallel_with;
        is_dependent_ = is_dependent;
        return;
    }

    /**
     * @brief Set the Joint Angle object
     * @param joint_angle [rad]
     * @todo 角度の変更に伴う同次変換行列の更新
     */
    void setJointAngle(const double joint_angle)
    {

        joint_angle_ = joint_angle * (does_reverse_ ? -1.0 : 1.0);
        z_link_rotate_ = get3DRotationMatrix(Eigen::Vector3d(0, 0, joint_angle_));
        return;
    }

    /**
     * @brief Adden joint angle to parallel link
     * @param add_joint_angle [rad]
     * @warning この関数はその関節の現在角度に対して角度を加算する.平行リンクなので同じ角度で符号が逆が加算される
     * @warning 加算は、idが大きい方のリンク要素が行う。小さい方が大きい方に加算してはいけない
     */
    void addenParallelJointAngle(const uint32_t &id_adden_with, const double add_joint_angle)
    {
        if (id_parallel_with_ != id_adden_with)
        {
            throw std::runtime_error("[addenParallelJointAngle()] id_parallel_with_ is not equal to id_adden_with !");
        }
        if (!is_dependent_)
        {
            // 従属側で無ければ加算しない
            return;
        }
        else
        {
            joint_angle_ -= add_joint_angle * (does_reverse_ ? -1.0 : 1.0);
            z_link_rotate_ = get3DRotationMatrix(Eigen::Vector3d(0, 0, joint_angle_));
            return;
        }
        return;
    }

    Eigen::Matrix4d getTransformMatrix()
    {
        Eigen::Matrix4d transform_matrix = Eigen::Matrix4d::Identity();
        transform_matrix.block<3, 1>(0, 3) = x_link_length_.transpose();
        transform_matrix *= x_link_rotate_;
        Eigen::Matrix4d tmp = Eigen::Matrix4d::Identity();
        tmp.block<3, 1>(0, 3) = z_link_length_.transpose();
        transform_matrix *= tmp;
        transform_matrix *= z_link_rotate_;
        return transform_matrix;
    }

    /**
     * @brief DH記法におけるzの回転角を更新する
     * @param angle 回転角度[rad]
     */
    void updateJointAngle(const double angle)
    {
        z_link_rotate_(0, 2) = angle;
        return;
    }

    /**
     * @brief 逆運動学を解く
     * @param target_position 目標位置(x,y,z) [mm]
     * @todo 逆運動学の実装 向きを決めて、sin,cosを解くのやつを実装する。
     */
    void solveInverseKinematicsGeometry(const Eigen::Vector3d &target_position)
    {
    }
};

/**
 * @brief アームのリンクそれぞれの位置の結果
 */
struct ResultPosition
{
    std::vector<float> x;
    std::vector<float> y;
    std::vector<float> z;
    auto getEndPositionVec()
    {
        return Eigen::Vector3d(x.back(), y.back(), z.back());
    }
    auto getPositionVec(const size_t index)
    {
        return Eigen::Vector3d(x[index], y[index], z[index]);
    }
    void print()
    {
        std::cout << "^^^ ResultPosition  ^^^\n";
        std::cout << "x: ";
        for (auto &data : x)
        {
            std::cout << data << ", ";
        }
        std::cout << "\n";
        std::cout << "y: ";
        for (auto &data : y)
        {
            std::cout << data << ", ";
        }
        std::cout << "\n";
        std::cout << "z: ";
        for (auto &data : z)
        {
            std::cout << data << ", ";
        }
        std::cout << "\n";
    }
    double EndX() const
    {
        return x.back();
    }
    
    double EndY() const
    {
        return y.back();
    }
    
    double EndZ() const
    {
        return z.back();
    }
};

struct RobotModel
{
    std::shared_ptr<std::vector<StraightChainRobotModel>> model_ptr_;
    std::vector<StraightChainRobotModel> model_;
    RobotModel() : model_ptr_(nullptr), model_()
    {
        // ロボットモデルの初期化
        createRobotModel();
    }
    void createRobotModel()
    {
        model_.clear();
        uint32_t id = 0;
        // ロボットの座標系との変換のために存在.ロボットは正面方向が+xだが、この変換が無いと-xが正面方向になってしまう
        model_.push_back(StraightChainRobotModel(id++, 0, Eigen::Vector3d::Zero(), 0, Eigen::Vector3d(0, 0, M_PI / 2)));
        model_.back().setImagenary(true);
        // 腰ヨー
        model_.push_back(StraightChainRobotModel(id++, 0, Eigen::Vector3d::Zero(), -L01_Y_WAIST_YAW, Eigen::Vector3d::Zero()));
        // 腰付け根ロール
        model_.push_back(StraightChainRobotModel(id++, 0, Eigen::Vector3d(M_PI / 2.0, 0, 0), 0, Eigen::Vector3d::Zero()));
        model_.back().setReverse(true); //ここ逆にしたらxv_refが楽
        // 腰独立への変換(Imagenary Link)
        model_.push_back(StraightChainRobotModel(id++, 0, Eigen::Vector3d(0, -M_PI / 2.0, 0), 0, Eigen::Vector3d(0, 0, -M_PI / 2.0)));
        model_.back().setImagenary(true);
        // 腰独立ピッチ
        model_.push_back(StraightChainRobotModel(id++, 0, Eigen::Vector3d(0, 0, 0), 0, Eigen::Vector3d(0, 0, 0)));
        model_.back().setParallelWith(id, true);
        // 膝上ピッチまで
        model_.push_back(StraightChainRobotModel(id++, L34_Z_WAIST_PITCH_TO_UPPER_KNEE, Eigen::Vector3d(0, 0, 0), 0, Eigen::Vector3d(0, 0, 0)));
        model_.back().setParallelWith(id - 2, false);
        model_.back().setReverse(true);
        // 膝上と膝下の間
        model_.push_back(StraightChainRobotModel(id++, L45_Z_UPPER_KNEE_TO_LOWER_KNEE, Eigen::Vector3d(0, 0, 0), 0, Eigen::Vector3d(0, 0, 0)));
        model_.back().setParallelWith(id + 1, false);
        // 膝下から足首まで
        model_.push_back(StraightChainRobotModel(id++, L56_Z_LOWER_KNEE_TO_FOOT_PITCH, Eigen::Vector3d(0, 0, 0), 0, Eigen::Vector3d(0, 0, 0)));
        model_.back().setImagenary(true);
        // 足首のピッチ回転
        model_.push_back(StraightChainRobotModel(id++, 0, Eigen::Vector3d(0, 0, 0), 0, Eigen::Vector3d(0, 0, 0)));
        model_.back().setParallelWith(id - 3, true);
        // model_.back().setReverse(true); ここ無くしたら右脚はxv_ref.dの符号変えずに行けた
        // 足首のロール回転
        model_.push_back(StraightChainRobotModel(id++, 0, Eigen::Vector3d(M_PI / 2, 0, 0), 0, Eigen::Vector3d(0, 0, 0)));
        // 足首から足裏まで
        model_.push_back(StraightChainRobotModel(id++, L67_Z_FOOT_PITCH_TO_FOOT_BOARD, Eigen::Vector3d(0, 0, 0), 0, Eigen::Vector3d(0, 0, 0)));
        model_.back().setImagenary(true);
    }
    std::shared_ptr<std::vector<StraightChainRobotModel>> getPtr()
    {
        return model_ptr_;
    }

    /**
     * @brief Set the All Joint Angle object
     * @param joint_angles 各関節角度　[rad]
     */
    void setAllJointAngle(const std::vector<double> joint_angles)
    {
        auto itr = joint_angles.begin();
        for (auto &motor : model_)
        {
            if (motor.is_imagenary_)
            {
                continue;
            }
            motor.setJointAngle(*itr);
            itr++;
        }
        for (auto &motor : model_)
        {
            if (motor.isParallel() && motor.isDependent())
            {
                motor.addenParallelJointAngle(motor.id_parallel_with_, model_[motor.id_parallel_with_].joint_angle_);
            }
        }
        return;
    }

    void printAllJointAnglesDeg()
    {
        for (auto &motor : model_)
        {
            if (motor.is_imagenary_)
            {
                continue;
            }
            std::cout << "motor id <" << motor.motor_id_ << "> angle <" << motor.joint_angle_ * 180.0 / M_PI << ">" << std::endl;
        }
        return;
    }

    ResultPosition calcForwardKinematics(bool print = false)
    {
        if (model_.size() == 0)
        {
            throw std::runtime_error("model_ is empty!");
        }
        ResultPosition result;
        result.x.emplace_back(0); // 原点
        result.y.emplace_back(0); // 原点
        result.z.emplace_back(0); // 原点
        Eigen::Matrix4d res = Eigen::Matrix4d::Identity();
        for (auto &motor : model_)
        {
            auto tf = motor.getTransformMatrix();
            res *= tf;
            result.x.emplace_back(res(0, 3));
            result.y.emplace_back(res(1, 3));
            result.z.emplace_back(res(2, 3));
            if (print)
            {
                std::cout << "*** motor id " << (uint16_t)(motor.motor_id_) << " ***" << std::endl;
                std::cout << tf << std::endl;
                std::cout << "--------------------------" << std::endl;
                std::cout << res << std::endl;
            }
        }
        if (print)
        {
            std::cout << "^^^ result ^^^" << std::endl;
            std::cout << res << std::endl;
        }
        return result;
    }
};

#endif // !ROBOT_MODEL_HPP_