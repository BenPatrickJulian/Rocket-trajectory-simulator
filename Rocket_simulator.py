# -*- coding: utf-8 -*-
"""
Created on Tue Jun  2 15:23:40 2026

@author: n14123bj
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Things changeable

dt = 0.001
Diameter_of_main_tube = 4.18*10**-2
empty_mass = 160*10**-3
fuel_mass = 30*10**-3
Total_impulse = 57.1
Name_of_motor_file = 'AeroTech_F35W.csv'
Graph = 'position'
percentages_shown = True
overlay_openrocket = True
Openrocket_file_name = 'openrocket_data.csv'

# Constants

mass_density_sea_level = 1.225
g = 9.80665
molar_mass_of_dry_air = 0.028964
R = 8.31446
Temperature = 287
drag_coefficient = 0.578
reference_area = np.pi * (Diameter_of_main_tube/2)**2
total_mass = empty_mass + fuel_mass
sea_level_standard_pressure = 101325
sea_level_standard_temperature = 288.15
Temperature_lapse_rate = 0.0065

# Set Things to 0

time = 0
acceleration = 0
velocity = 0 
position = 0
time_array = []
position_array = []
velocity_array = []
acceleration_array = []
position_comparison_array = [0,0,0]
velocity_comparison_array = [0,0,0]
Impulse_at_time_t = 0
F_previous = 0
F_thrust = 0
max_velocity = 0
max_acceleration = 0
apogee = 0
apogee_time = 0
m_air_density = 0
Graph = Graph.lower()

# Functions

def get_data_from_file(file_name):
    try:
        raw_data = pd.read_csv(file_name, comment='#').to_numpy().astype(object)
    except FileNotFoundError:
        print("File with given name does not exist in the directory")
        sys.exit()
    clean_rows = []

    for row in raw_data:
        redact_line = False
        converted_row = []

        for item in row:
            if pd.isna(item):
                redact_line = True
            try:
                value = float(item)
                converted_row.append(value)
            except (ValueError, TypeError):
                redact_line = True

        if not redact_line:
            clean_rows.append(np.array(converted_row))

    return np.array(clean_rows)

def calculate_drag(mass_density, drag_c, ref_area, velocity):
    F_drag = 0.5*mass_density*drag_c*ref_area*(velocity**2)*np.sign(velocity)
    return F_drag

def update_mass(empty_mass, fuel_mass, time, Area_under_curve, Total_impulse):
    total_mass = empty_mass + fuel_mass - (Area_under_curve/Total_impulse)*fuel_mass
    return total_mass

def update_acceleration(F_drag, F_thrust, total_mass, fuel_empty):
    F_weight = total_mass*9.81
    if fuel_empty == True:
        F_thrust = 0
    F_total = F_thrust - F_drag - F_weight
    acceleration = F_total/total_mass
    return acceleration

def update_velocity(acceleration, dt, velocity):
    velocity = velocity + acceleration*dt
    return velocity

def update_position(acceleration, velocity, dt, position):
    position = position + velocity*dt + 0.5*acceleration*(dt**2)
    return position

def check_fuel_supply(Total_impulse, Area_under_curve):
    if Area_under_curve/Total_impulse >= 0.99:
        return True
    else:
        return False
    
def update_apogee(position, apogee, time, apogee_time):
    if position >= apogee:
        apogee = position
        apogee_time = time
    return apogee, apogee_time

def update_max_velocity(velocity, max_velocity):
    if velocity >= max_velocity:
        max_velocity = velocity
    return max_velocity

def update_max_acceleration(acceleration, max_acceleration):
    if acceleration >= max_acceleration:
        max_acceleration = acceleration
    return max_acceleration

def update_Thrust_and_impulse(F_thrust, time, time_values, thrust_values, Impulse_at_time_t):
    F_previous = F_thrust
    F_thrust = np.interp(time, time_values, thrust_values, right=0.0)
    Impulse_added = dt*0.5*(F_thrust + F_previous)
    Impulse_at_time_t = Impulse_at_time_t + Impulse_added
    return F_thrust, Impulse_at_time_t

def update_mass_density(sea_level_standard_pressure, Temperature_lapse_rate, sea_level_standard_temperature, position):
    pressure = sea_level_standard_pressure*(1-(Temperature_lapse_rate*position)/sea_level_standard_temperature)**(5.25588)
    mass_density = pressure/(287.05*(sea_level_standard_temperature-Temperature_lapse_rate*position))
    return mass_density

def calculate_NRMSE(time_array, position_array, velocity_array, acceleration_array, OR_position_array, OR_velocity_array, OR_acceleration_array, apogee, max_velocity, max_acceleration):
    position_sum = 0
    velocity_sum = 0
    acceleration_sum = 0
    for i in range(len(time_array)):
        position_value = (position_array[i] - OR_position_array[i])**2
        position_sum = position_sum + position_value
        velocity_value = (velocity_array[i] - OR_velocity_array[i])**2
        velocity_sum = velocity_sum + velocity_value
        acceleration_value = (acceleration_array[i] - OR_acceleration_array[i])**2
        acceleration_sum = acceleration_sum + acceleration_value
        
    NRMSE_position = np.sqrt(1/len(time_array) * position_sum)/apogee
    NRMSE_velocity = np.sqrt(1/len(time_array) * velocity_sum)/max_velocity
    NRMSE_acceleration = np.sqrt(1/len(time_array) * acceleration_sum)/max_acceleration
    
    print("The position was", 100 - NRMSE_position*100, "% correct")
    print("The velocity was", 100 - NRMSE_velocity*100, "% correct")
    print("The acceleration was", 100 - NRMSE_acceleration*100, "% correct")

def show_openrocket_overlay(time_array, position_array, velocity_array, OR_position_array, OR_velocity_array):
    if Graph == 'position':
        plt.plot(time_array, position_array, label='Data')
        plt.plot(time_array, OR_position_array, label='Openrocket Data')
        plt.ylabel("Height/m")
    elif Graph == 'velocity':
        plt.plot(time, velocity_array, label='Data')
        plt.plot(time, OR_velocity_array, label='Openrocket Data')
        plt.ylabel("Velocity/ms-1")
    plt.xlabel("Time/s")
    plt.legend()
    plt.grid()
 
thrust_array = get_data_from_file(Name_of_motor_file)
time_values = thrust_array[:, 0]
thrust_values = thrust_array[:, 1]

Flying = True
while Flying == True:
    time = time + dt
    fuel_empty = check_fuel_supply(Total_impulse, Impulse_at_time_t)
    F_thrust, Impulse_at_time_t = update_Thrust_and_impulse(F_thrust, time, time_values, thrust_values, Impulse_at_time_t)
    if fuel_empty == False:
        total_mass = update_mass(empty_mass, fuel_mass, time, Impulse_at_time_t, Total_impulse)
    air_mass_density = update_mass_density(sea_level_standard_pressure, Temperature_lapse_rate, sea_level_standard_temperature, position)
    if position >= 8000 and position <= 8020:
        m_air_density = air_mass_density
    F_drag = calculate_drag(air_mass_density, drag_coefficient, reference_area, velocity)
    acceleration = update_acceleration(F_drag, F_thrust, total_mass, fuel_empty)
    velocity = update_velocity(acceleration, dt, velocity)
    position = update_position(acceleration, velocity, dt, position)
    if position<-0.01:
        Flying = False 
    apogee, apogee_time = update_apogee(position, apogee, time, apogee_time)
    max_velocity = update_max_velocity(velocity, max_velocity)
    max_acceleration = update_max_acceleration(acceleration, max_acceleration)
    time_array.append(time)
    position_array.append(position)
    velocity_array.append(velocity)
    acceleration_array.append(acceleration)

print(f"Apogee was {apogee:.3g}m at {apogee_time:.3g}s")
print(f"Max velocity was {max_velocity:.3g}ms-1")
print(f"Max acceleration was {max_acceleration:.3g}ms-1")

if percentages_shown == True:
    openrocket_data = get_data_from_file(Openrocket_file_name)
    accurate_time = openrocket_data[:, 0]
    accurate_position = openrocket_data[:, 1]
    accurate_velocity = openrocket_data[:, 2]
    accurate_acceleration = openrocket_data[:, 3]
    interp_accurate_position = np.interp(time_array, accurate_time, accurate_position)
    interp_accurate_velocity = np.interp(time_array, accurate_time, accurate_velocity)
    interp_accurate_acceleration = np.interp(time_array, accurate_time, accurate_acceleration)
    calculate_NRMSE(time_array, position_array, velocity_array, acceleration_array, interp_accurate_position, interp_accurate_velocity, interp_accurate_acceleration, apogee, max_velocity, max_acceleration)

if overlay_openrocket == True:
    show_openrocket_overlay(time_array, position_array, velocity_array, interp_accurate_position, interp_accurate_velocity)

if Graph == 'position':
    plt.plot(time_array, position_array)
    plt.ylabel("Height/m")
elif Graph == 'velocity':
    plt.plot(time_array, velocity_array)
    plt.ylabel("Velocity/ms-1")
plt.xlabel("Time/s")
plt.grid()
plt.show()    
    





