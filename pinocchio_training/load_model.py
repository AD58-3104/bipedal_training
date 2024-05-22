import pinocchio
from sys import argv
from os.path import dirname, join, abspath
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# You should change here to set up your own URDF file or just pass it as an argument of this example.
urdf_filename = argv[1]

# Load the urdf model
model    = pinocchio.buildModelFromUrdf(urdf_filename)
print('model name: ' + model.name)

# Create data required by the algorithms
data     = model.createData()

# Sample a random configuration
q        = pinocchio.randomConfiguration(model)
print('q: %s' % q.T)

# Perform the forward kinematics over the kinematic tree
pinocchio.forwardKinematics(model,data,q)

# Print out the placement of each joint of the kinematic tree
for name, oMi in zip(model.names, data.oMi):
    print(("{:<24} : {: .2f} {: .2f} {: .2f}"
          .format( name, *oMi.translation.T.flat )))


fig = plt.figure()
# Create a figure and a 3D axis
ax = fig.add_subplot(111, projection='3d')
# Plot the arm joints
for i in range(len(model.names)-1):
    x1, y1, z1 = data.oMi[i].translation.T.flat
    x2, y2, z2 = data.oMi[i+1].translation.T.flat
    ax.plot([x1, x2], [y1, y2], [z1, z2], label=model.names[i])
    # Plot the last point as a larger, red sphere
    x_last, y_last, z_last = data.oMi[-1].translation.T.flat
ax.scatter(x_last, y_last, z_last, color='red', s=100)
# Set axis labels
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

# Set plot title
ax.set_title('Arm Position')

# Show the plot
plt.show()