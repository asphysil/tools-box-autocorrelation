import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

# Read x, y dataset from file
with open('ref-test-vac.dat', 'r') as file:
    data = file.readlines()

x_data = []
y_data = []

for line in data:
    values = line.strip().split()  # Split by any whitespace (spaces, tabs, etc.)
    x, y = map(float, values)
    x_data.append(x)
    y_data.append(y)

# Convert lists to NumPy arrays
x = np.array(x_data)
y = np.array(y_data)

# Interpolate to create more data points
f = interp1d(x, y, kind='linear')  # You can choose 'linear', 'quadratic', 'cubic', etc.
new_x = np.linspace(min(x), max(x), num=6000)  # Increase the number of data points
new_y = f(new_x)

tmax = max(new_x) 

tdamp = tmax/np.sqrt(3.0)

new_y2 = new_y*np.exp(-(new_x/tdamp)**2) 

m = 2**(int(np.log(len(new_x))/np.log(2.0)))

print('number of data set', m, len(new_x))

n = len(new_x)
 
w = np.fft.fft(new_y2 - np.mean(new_y2), n=2*n)

w_norm = w/n 

freq = np.fft.fftfreq(2*n, d=new_x[1] - new_x[0])
magnitude_spectrum = np.abs(w_norm)**2

print(freq[:int(len(freq)/2)], magnitude_spectrum)

plt.figure(figsize=(8, 6))
plt.plot(freq[:int(len(freq)/2)]*2*np.pi, magnitude_spectrum[:int(len(freq)/2)])
plt.title('Magnitude Spectrum')
plt.xlabel('Frequency')
plt.ylabel('Magnitude')
plt.grid(True)
plt.xlim([0,  30])


# Plotting the original and interpolated data
plt.figure(figsize=(8, 6))
plt.scatter(x, y, color='blue', label='Original Data')
plt.plot(new_x, new_y, color='red', label='Interpolated Data')
plt.xlabel('X-axis')
plt.ylabel('Y-axis')
plt.title('Original and Interpolated Data')
plt.legend()
plt.grid(True)
plt.show()

# Write the expanded dataset to a new file
with open('expanded_data.txt', 'w') as output_file:
    for i in range(len(new_x)):
        output_file.write(f"{new_x[i]}  {new_y[i]}\n")

print("Data has been written to expanded_data.txt")

