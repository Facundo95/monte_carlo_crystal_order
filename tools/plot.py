#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 19 17:27:03 2026

@author: facundo
"""
import matplotlib.pyplot as plt
import numpy as np
import argparse

def plot(temp, campo, lro, magn, etotal):

    fig, ax = plt.subplots()
    for t, p in zip(temp, lro):
        ax.plot(t, p)
    ax.set_xlabel('Temperature (K)')
    ax.set_ylabel('lro Parameters (u.a.)')
    ax.legend()

    ax2 = ax.twinx()
    ax2.plot(temp, magn)
    ax2.set_ylabel('Magnetization (u.a.)')
    plt.show()

def get_data(file_path):
    return np.loadtxt(file_path, header=1, delimiter='\t')

def process_data(data):
    temp = data[:, 0]
    campo = data[:, 1]
    lro = data[:, 2:5]  # Assuming lro parameters are in columns 2, 3, and 4
    magn = data[:, 5]    # Assuming magnetization is in column 5
    etotal = data[:, 6] # Assuming total energy is in column 6
    return temp, campo, lro, magn, etotal

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot lro parameters and magnetization.')
    parser.add_argument('-f', type=str, help='Path to the data file .out containing temperature, lro parameters, and magnetization.')
    args = parser.parse_args()

    data = get_data(args.f)
    temp, campo, lro, magn, etotal = process_data(data)
    plot(temp, campo, lro, magn, etotal)