import zmq

import numpy as np
import matplotlib.pyplot as plt

ctx = zmq.Context()
socket = ctx.socket(zmq.SUB)
socket.connect( "tcp://127.0.0.1:9000")
print("SUB connected")
#socket.setsockopt(zmq.SUBSCRIBE, b"")
fig, ax = plt.subplots(2, 2, figsize=(10, 6))
ax1 = ax[0, 0]
ax2 = ax[0, 1]
ax3 = ax[1, 0]
ax4 = ax[1, 1]
L = 30
rc = 2.5
ax1.grid(True, which='both', linestyle='-', linewidth=1)
ax1.set_xticks(np.arange(0, L + rc, rc))
ax1.set_yticks(np.arange(0, L + rc, rc))
ax1.set_xlim([0, L])
ax1.set_ylim([0, L])
potential_energy_list = []
kinetic_energy_list = []
total_energy_list = []
momentum_list = []
temperature_list = []
scatter_plot = ax1.scatter([], [], c='r')
#pot_en_line, = ax2.plot([], [], c='red')
#kin_en_line, = ax2.plot([], [], c='blue')
#tot_en_line, = ax2.plot([], [], c='green')
#mom_line, = ax3.plot([], [], c='green')
#temp_line, = ax4.plot([], [], c='red')
#socket.setsockopt( zmq.LINGER, 0 )
#socket.setsockopt( zmq.SUBSCRIBE, "" )

line0 = ax2.plot(potential_energy_list, c='blue')
line1 = ax2.plot(kinetic_energy_list, c='red')
line2 = ax2.plot(total_energy_list, c='green')
line3 = ax3.plot(momentum_list, c='red')
line4 = ax4.plot(temperature_list, c='blue')
print("start receiving")
while 1:


    #frames = socket.recv_multipart()  # Receive multiple frames
    #print(frames)  # Print all frames in the message
    x = np.frombuffer(socket.recv(), dtype=np.float64)
    print(x)
    y = np.frombuffer(socket.recv(), dtype=np.float64)
    pot_en = np.frombuffer(socket.recv(), dtype=np.float64)

    kin_en = np.frombuffer(socket.recv(), dtype=np.float64)
    momentum = np.frombuffer(socket.recv(), dtype=np.float64)
    T = np.frombuffer(socket.recv(), dtype=np.float64)
    #print("Received request: ", x, y, pot_en, kin_en, momentum, T)
    potential_energy_list.append(pot_en)
    kinetic_energy_list.append(kin_en)
    total_energy_list.append(pot_en + kin_en)
    momentum_list.append(momentum)
    temperature_list.append(T)
    scatter_plot.set_offsets(np.column_stack((x, y)))
    print("messages received")
    #socket.send(b"OK")
    line0.set_ydata(kinetic_energy_list)
    line1.set_ydata(potential_energy_list)
    line2.set_ydata(momentum_list)
    line3.set_ydata(temperature_list)

    ax2.set_xlabel('Iter')
    ax2.set_ylabel('Energy')
    ax3.set_xlabel('Iter')
    ax3.set_ylabel('Momentum')
    ax4.set_xlabel('Iter')
    ax4.set_ylabel('Temperature')
    plt.draw()
    plt.pause(0.000001)
socket.disconnect()




    #pot_en_line.set_data(list(range(len(potential_energy_list))), potential_energy_list)
    #kin_en_line.set_data(list(range(len(kinetic_energy_list))), kinetic_energy_list)
    #tot_en_line.set_data(list(range(len(total_energy_list))), total_energy_list)
    #mom_line.set_data(list(range(len(momentum_list))), momentum_list)
    #temp_line.set_data(list(range(len(temperature_list))), temperature_list)

    #fig.canvas.draw()


