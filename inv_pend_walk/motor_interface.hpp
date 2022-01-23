#ifndef WEBOTS_INTERFACE
#define WEBOTS_INTERFACE
#include <webots/Motor.hpp>
#include <webots/Robot.hpp>
#include <webots/Accelerometer.hpp>
#include <webots/Gyro.hpp>
#include <unordered_map>
#include <set>
enum class motor_label
{
    FOOT_ROLL_R,
    LEG_PITCH_R,
    KNEE_R1,
    KNEE_R2,
    LEG_ROLL_R,
    LEG_YAW_R,
    ARM_ROLL_R,
    ARM_PITCH_R,
    ELBOW_PITCH_R,
    FOOT_ROLL_L,
    LEG_PITCH_L,
    KNEE_L1,
    KNEE_L2,
    LEG_ROLL_L,
    LEG_YAW_L,
    ARM_PITCH_L,
    ARM_ROLL_L,
    ELBOW_PITCH_L,
    HEAD_YAW
};

class webots_control
{
private:
    std::unordered_map<motor_label, std::string> motor_names; // motor_names
    webots::Robot *robot;
    std::unordered_map<motor_label, std::pair<webots::Motor *, std::string>> robot_motors;
    std::set<std::string> reverse_motors;
    int32_t mTimeStep;
public:
    webots_control() : mTimeStep(0)
    {
        robot = new webots::Robot();
        motor_names[motor_label::FOOT_ROLL_R] =  "right_ankle_roll_joint";
        motor_names[motor_label::LEG_PITCH_R] =  "right_ankle_pitch_joint";
        motor_names[motor_label::KNEE_R1] =  "right_knee_pitch_joint";
        motor_names[motor_label::KNEE_R2] =  "right_waist_pitch_joint";
        motor_names[motor_label::LEG_ROLL_R] =  "right_waist_roll_joint [hip]";
        motor_names[motor_label::LEG_YAW_R] =  "right_waist_yaw_joint";
        motor_names[motor_label::ARM_ROLL_R] =  "right_shoulder_roll_joint";
        motor_names[motor_label::ARM_PITCH_R] =  "right_shoulder_pitch_joint [shoulder]";
        motor_names[motor_label::ELBOW_PITCH_R] =  "right_elbow_pitch_joint";
        motor_names[motor_label::FOOT_ROLL_L] =  "left_ankle_roll_joint";
        motor_names[motor_label::LEG_PITCH_L] =  "left_ankle_pitch_joint";
        motor_names[motor_label::KNEE_L1] =  "left_knee_pitch_joint";
        motor_names[motor_label::KNEE_L2] =  "left_waist_pitch_joint";
        motor_names[motor_label::LEG_ROLL_L] =  "left_waist_roll_joint [hip]";
        motor_names[motor_label::LEG_YAW_L] =  "left_waist_yaw_joint";
        motor_names[motor_label::ARM_PITCH_L] =  "left_shoulder_pitch_joint [shoulder]";
        motor_names[motor_label::ARM_ROLL_L] =  "left_shoulder_roll_joint";
        motor_names[motor_label::ELBOW_PITCH_L] =  "left_elbow_pitch_joint";
        motor_names[motor_label::HEAD_YAW] =  "head_yaw_joint";

        // reverse_motors.emplace("right_knee_pitch_mimic_joint");
        // reverse_motors.emplace("right_ankle_pitch_mimic_joint");
        // reverse_motors.emplace("left_knee_pitch_mimic_joint");
        // reverse_motors.emplace("left_ankle_pitch_mimic_joint");
        reverse_motors.emplace("right_knee_pitch_joint");
        reverse_motors.emplace("right_waist_pitch_joint");
        reverse_motors.emplace("left_knee_pitch_joint");
        reverse_motors.emplace("left_waist_pitch_joint");
        reverse_motors.emplace("right_shoulder_roll_joint");
        reverse_motors.emplace("left_shoulder_roll_joint");
    }
};

#endif // !WEBOTS_INTERFACE