import zmq
import numpy as np
import matplotlib.pyplot as plt

ctx = zmq.Context()
socket = ctx.socket(zmq.SUB)
socket.connect("tcp://localhost:5556")
socket.setsockopt(zmq.SUBSCRIBE, b"")
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


while 1:
    x = np.frombuffer(socket.recv(), dtype=np.float64)
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
    #socket.send(b"OK")
    ax2.plot(list(range(len(potential_energy_list))), potential_energy_list, c='blue')
    ax2.plot(list(range(len(kinetic_energy_list))), kinetic_energy_list, c='red')
    ax2.plot(list(range(len(total_energy_list))), total_energy_list, c='green')
    ax3.plot(list(range(len(momentum_list))), momentum_list, c='red')
    ax4.plot(list(range(len(temperature_list))), temperature_list, c='blue')
    ax2.set_xlabel('Iter')
    ax2.set_ylabel('Energy')
    ax3.set_xlabel('Iter')
    ax3.set_ylabel('Momentum')
    ax4.set_xlabel('Iter')
    ax4.set_ylabel('Temperature')
    plt.pause(0.0001)





    #pot_en_line.set_data(list(range(len(potential_energy_list))), potential_energy_list)
    #kin_en_line.set_data(list(range(len(kinetic_energy_list))), kinetic_energy_list)
    #tot_en_line.set_data(list(range(len(total_energy_list))), total_energy_list)
    #mom_line.set_data(list(range(len(momentum_list))), momentum_list)
    #temp_line.set_data(list(range(len(temperature_list))), temperature_list)

    #fig.canvas.draw()


