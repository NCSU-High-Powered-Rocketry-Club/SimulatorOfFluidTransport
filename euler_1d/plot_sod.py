import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("sod_results.csv")

plt.figure(figsize=(10, 6))
plt.plot(df['x'], df['density'], label='Density', color='black')
plt.plot(df['x'], df['pressure'], label='Pressure', linestyle='--', color='red')
plt.plot(df['x'], df['velocity'], label='Velocity', linestyle=':', color='blue')
plt.title("Sod Shock Tube: 1st Order Euler with Edwards LDFSS")
plt.legend()
plt.grid(True)
plt.show()
