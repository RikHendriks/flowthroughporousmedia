from core import *


# Instantiate an instance of the map finite difference class

f = MapFinDif(5.0, 5.0, 0.1)


# Setting boundary conditions

# Horizontal plane y=0
for c in f.horizontal_up_cell_generator([0.0, 5.0], 0.0):
    c.bound_down = 0.0

# Horizontal plane y=3
for c in f.horizontal_down_cell_generator([0.0, 5.0], 5.0):
    c.bound_up = 3.0

# Vertical plane x=0
for c in f.vertical_right_cell_generator([0.0, 5.0], 0.0):
    c.bound_left = 1.0

# Vertical plane x=20
for c in f.vertical_left_cell_generator([0.0, 5.0], 5.0):
    c.bound_right = 2.0


# Calculations

f.optimize(0.5, 0.0001)


# Visualization

f.show_discharge()
