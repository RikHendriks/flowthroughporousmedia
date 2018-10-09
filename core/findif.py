import numpy as np
from matplotlib import pyplot as plt


class MapFinDif:
    def __init__(self, width, height, element_size=0.1, hydraulic_conductivity=1.0):
        # Map properties
        self.width = width
        self.height = height
        # Element properties
        self.element_size = element_size
        self.x_elements = int(self.width / self.element_size)
        self.y_elements = int(self.height / self.element_size)
        self.map = [[Cell([i, j], np.random.random_sample()) for j in range(self.y_elements)] for i in range(self.x_elements)]
        # Properties
        self.hydraulic_conductivity = hydraulic_conductivity

    def laplace_finite_difference_cell(self, coordinate):
        x_neighbour_value = 0.0
        y_neighbour_value = 0.0
        laplace_factor = 4.0

        # Neumann boundary conditions x axis
        if self.map[coordinate[0]][coordinate[1]].bound_flow_left:
            laplace_factor -= 1.0
        # Dirichlet boundary conditions x axis left
        elif self.map[coordinate[0]][coordinate[1]].bound_left is not None:
            x_neighbour_value += self.map[coordinate[0]][coordinate[1]].bound_left
        else:
            x_neighbour_value += self.map[coordinate[0] - 1][coordinate[1]].value

        # Neumann boundary conditions x axis
        if self.map[coordinate[0]][coordinate[1]].bound_flow_right:
            laplace_factor -= 1.0
        # Dirichlet boundary conditions x axis right
        elif self.map[coordinate[0]][coordinate[1]].bound_right is not None:
            x_neighbour_value += self.map[coordinate[0]][coordinate[1]].bound_right
        else:
            x_neighbour_value += self.map[coordinate[0] + 1][coordinate[1]].value

        # Neumann boundary conditions x axis
        if self.map[coordinate[0]][coordinate[1]].bound_flow_down:
            laplace_factor -= 1.0
        # Dirichlet Boundary conditions y axis down
        elif self.map[coordinate[0]][coordinate[1]].bound_down is not None:
            y_neighbour_value += self.map[coordinate[0]][coordinate[1]].bound_down
        else:
            y_neighbour_value += self.map[coordinate[0]][coordinate[1] - 1].value

        # Neumann boundary conditions x axis
        if self.map[coordinate[0]][coordinate[1]].bound_flow_up:
            laplace_factor -= 1.0
        # Dirichlet Boundary conditions y axis up
        elif self.map[coordinate[0]][coordinate[1]].bound_up is not None:
            y_neighbour_value += self.map[coordinate[0]][coordinate[1]].bound_up
        else:
            y_neighbour_value += self.map[coordinate[0]][coordinate[1] + 1].value

        # Return the value of the laplace finite difference
        if laplace_factor == 0.0:
            return 0.0
        else:
            return (1 / laplace_factor) * (x_neighbour_value + y_neighbour_value)

    def derivative_x_cell(self, coordinate):
        x_0_neighbour_value = 0.0
        x_1_neighbour_value = 0.0
        first = True
        second = True

        # Neumann boundary conditions x axis
        if self.map[coordinate[0]][coordinate[1]].bound_flow_left:
            first = False
        # Dirichlet boundary conditions x axis left
        elif self.map[coordinate[0]][coordinate[1]].bound_left is not None:
            x_0_neighbour_value += self.map[coordinate[0]][coordinate[1]].bound_left
        else:
            x_0_neighbour_value += self.map[coordinate[0] - 1][coordinate[1]].value

        # Neumann boundary conditions x axis
        if self.map[coordinate[0]][coordinate[1]].bound_flow_right:
            second = False
        # Dirichlet boundary conditions x axis right
        elif self.map[coordinate[0]][coordinate[1]].bound_right is not None:
            x_1_neighbour_value += self.map[coordinate[0]][coordinate[1]].bound_right
        else:
            x_1_neighbour_value += self.map[coordinate[0] + 1][coordinate[1]].value

        # Central difference
        if first and second:
            return (x_1_neighbour_value - x_0_neighbour_value) / (2 * self.element_size)
        # Backward difference
        elif first:
            return (self.map[coordinate[0]][coordinate[1]].value - x_0_neighbour_value) / self.element_size
        # Forward difference
        elif second:
            return (x_1_neighbour_value - self.map[coordinate[0]][coordinate[1]].value) / self.element_size

    def derivative_y_cell(self, coordinate):
        y_0_neighbour_value = 0.0
        y_1_neighbour_value = 0.0
        first = True
        second = True

        # Neumann boundary conditions x axis
        if self.map[coordinate[0]][coordinate[1]].bound_flow_down:
            first = False
        # Dirichlet Boundary conditions y axis down
        elif self.map[coordinate[0]][coordinate[1]].bound_down is not None:
            y_0_neighbour_value += self.map[coordinate[0]][coordinate[1]].bound_down
        else:
            y_0_neighbour_value += self.map[coordinate[0]][coordinate[1] - 1].value

        # Neumann boundary conditions x axis
        if self.map[coordinate[0]][coordinate[1]].bound_flow_up:
            second = False
        # Dirichlet Boundary conditions y axis up
        elif self.map[coordinate[0]][coordinate[1]].bound_up is not None:
            y_1_neighbour_value += self.map[coordinate[0]][coordinate[1]].bound_up
        else:
            y_1_neighbour_value += self.map[coordinate[0]][coordinate[1] + 1].value

        # Central difference
        if first and second:
            return (y_1_neighbour_value - y_0_neighbour_value) / (2 * self.element_size)
        # Backward difference
        elif first:
            return (self.map[coordinate[0]][coordinate[1]].value - y_0_neighbour_value) / self.element_size
        # Forward difference
        elif second:
            return (y_1_neighbour_value - self.map[coordinate[0]][coordinate[1]].value) / self.element_size

    def calculate_total_loss(self):
        total_loss = 0.0
        for i in range(self.x_elements):
            for j in range(self.y_elements):
                if self.map[i][j].active:
                    total_loss += (self.map[i][j].value - self.laplace_finite_difference_cell([i, j])) ** 2.0
        # Return the total loss
        return total_loss

    def coordinate_array(self):
        x_array = [(i + 0.5) * self.element_size for i in range(self.x_elements)]
        y_array = [(i + 0.5) * self.element_size for i in range(self.y_elements)]
        # Return the coordinate array
        return x_array, y_array

    def get_value_array(self):
        v_array = []
        for j in range(self.y_elements):
            v_array.append([])
            for i in range(self.x_elements):
                if self.map[i][j].active:
                    v_array[j].append(self.map[i][j].value)
                else:
                    v_array[j].append(0.0)
        # Return the v_array
        return v_array

    def get_discharge_x(self):
        d_x_array = []
        for j in range(self.y_elements):
            d_x_array.append([])
            for i in range(self.x_elements):
                if self.map[i][j].active:
                    d_x = -self.hydraulic_conductivity * self.derivative_x_cell([i, j])
                else:
                    d_x = 0.0
                d_x_array[j].append(d_x)
        # Return the d_x_array
        return d_x_array

    def get_discharge_y(self):
        d_y_array = []
        for j in range(self.y_elements):
            d_y_array.append([])
            for i in range(self.x_elements):
                if self.map[i][j].active:
                    d_y = -self.hydraulic_conductivity * self.derivative_y_cell([i, j])
                else:
                    d_y = 0.0
                d_y_array[j].append(d_y)
        # Return the d_y_array
        return d_y_array

    def show_head(self):
        x, y = self.coordinate_array()
        z = self.get_value_array()

        plt.contourf(x, y, z, 25)
        plt.colorbar()
        plt.show()

    def show_discharge_x(self):
        x, y = self.coordinate_array()
        z = self.get_discharge_x()

        plt.contourf(x, y, z, 25)
        plt.colorbar()
        plt.show()

    def show_discharge_y(self):
        x, y = self.coordinate_array()
        z = self.get_discharge_y()

        plt.contourf(x, y, z, 25)
        plt.colorbar()
        plt.show()

    def show_discharge(self):
        x, y = self.coordinate_array()
        u = self.get_discharge_x()
        v = self.get_discharge_y()

        plt.quiver(x, y, u, v)
        plt.show()

    def update_cells(self, over_relaxation):
        for i in range(self.x_elements):
            for j in range(self.y_elements):
                if self.map[i][j].active:
                    self.map[i][j].value += (over_relaxation + 1) * \
                                            (self.laplace_finite_difference_cell([i, j]) - self.map[i][j].value)

    def optimize(self, over_relaxation, error=0.001):
        while self.calculate_total_loss() > error:
            self.update_cells(over_relaxation)

    def square_cell_generator(self, coordinate_0, coordinate_1):
        for i in range(self.x_elements):
            for j in range(self.y_elements):
                coordinate = [(i + 0.5) * self.element_size,
                              (j + 0.5) * self.element_size]
                # If the coordinate is within the square then yield it
                if (coordinate_0[0] < coordinate[0] < coordinate_1[0] or
                    coordinate_0[0] > coordinate[0] > coordinate_1[0]) and \
                        (coordinate_0[1] < coordinate[1] < coordinate_1[1] or
                         coordinate_0[1] > coordinate[1] > coordinate_1[1]):
                    yield self.map[i][j]

    def horizontal_down_cell_generator(self, x_coordinates, y_axis):
        for j in range(self.y_elements):
            y = (j + 0.5) * self.element_size
            if y + self.element_size >= y_axis:
                for i in range(self.x_elements):
                    x = (i + 0.5) * self.element_size
                    if x_coordinates[0] <= x <= x_coordinates[1] or \
                            x_coordinates[0] >= x >= x_coordinates[1]:
                        yield self.map[i][j]
                break

    def horizontal_up_cell_generator(self, x_coordinates, y_axis):
        for j in range(self.y_elements):
            y = (j + 0.5) * self.element_size
            if y >= y_axis:
                for i in range(self.x_elements):
                    x = (i + 0.5) * self.element_size
                    if x_coordinates[0] <= x <= x_coordinates[1] or \
                            x_coordinates[0] >= x >= x_coordinates[1]:
                        yield self.map[i][j]
                break

    def vertical_left_cell_generator(self, y_coordinates, x_axis):
        for i in range(self.x_elements):
            x = (i + 0.5) * self.element_size
            if x + self.element_size >= x_axis:
                for j in range(self.y_elements):
                    y = (j + 0.5) * self.element_size
                    if y_coordinates[0] <= y <= y_coordinates[1] or \
                            y_coordinates[0] >= y >= y_coordinates[1]:
                        yield self.map[i][j]
                break

    def vertical_right_cell_generator(self, y_coordinates, x_axis):
        for i in range(self.x_elements):
            x = (i + 0.5) * self.element_size
            if x >= x_axis:
                for j in range(self.y_elements):
                    y = (j + 0.5) * self.element_size
                    if y_coordinates[0] <= y <= y_coordinates[1] or \
                            y_coordinates[0] >= y >= y_coordinates[1]:
                        yield self.map[i][j]
                break


class Cell:
    def __init__(self, coordinate, value,
                 bound_up=None, bound_down=None, bound_right=None, bound_left=None,
                 bound_flow_up=False, bound_flow_down=False, bound_flow_right=False, bound_flow_left=False):
        self.coordinate = coordinate
        self.value = value
        # Active
        self.active = True
        # Dirichlet boundaries
        self.bound_up = bound_up
        self.bound_down = bound_down
        self.bound_right = bound_right
        self.bound_left = bound_left
        # Neumann boundaries
        self.bound_flow_up = bound_flow_up
        self.bound_flow_down = bound_flow_down
        self.bound_flow_right = bound_flow_right
        self.bound_flow_left = bound_flow_left
