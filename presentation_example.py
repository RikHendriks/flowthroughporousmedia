from core import *


# Instantiate an instance of the map finite difference class

f = MapFinDif(20.0, 8.0, 0.25, 8e-5)


# Setting boundary conditions

# Inactive area
for c in f.square_cell_generator([10.0, 3.0], [15.0, 8.0]):
    c.active = False

# Horizontal plane y=0
for c in f.horizontal_up_cell_generator([0.0, 20.0], 0.0):
    c.bound_flow_down = True

# Horizontal plane y=3
for c in f.horizontal_down_cell_generator([10.0, 15.0], 3.0):
    c.bound_up = 3.0

# Horizontal plane y=8
for c in f.horizontal_down_cell_generator([0.0, 10.0], 8.0):
    c.bound_up = 8.0
for c in f.horizontal_down_cell_generator([15.0, 20.0], 8.0):
    c.bound_up = 8.0

# Vertical plane x=0
for c in f.vertical_right_cell_generator([0.0, 8.0], 0.0):
    c.bound_left = 8.0

# Vertical plane x=10
for c in f.vertical_left_cell_generator([2.0, 8.0], 10.0):
    c.bound_flow_right = True
for c in f.vertical_right_cell_generator([2.0, 8.0], 10.0):
    c.bound_flow_left = True

# Vertical plane x=15
for c in f.vertical_left_cell_generator([2.0, 8.0], 15.0):
    c.bound_flow_right = True
for c in f.vertical_right_cell_generator([2.0, 8.0], 15.0):
    c.bound_flow_left = True

# Vertical plane x=20
for c in f.vertical_left_cell_generator([0.0, 8.0], 20.0):
    c.bound_flow_right = True


# Calculations

f.optimize(0.5, 0.00001)


# Visualization

f.show_discharge()
