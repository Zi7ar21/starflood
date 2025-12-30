import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

log_dir = './out'

#log_name = 'log_nvidia-smi'
log_name = 'statistics'
#log_name = 'timings_sim_step'
#log_name = 'timings_vis_draw'

log_file = log_dir + '/' + log_name + '.csv'

if ('log_nvidia-smi' == log_name):
	df = pd.read_csv(log_file, parse_dates=['timestamp'])
else:
	df = pd.read_csv(log_file)

df.columns = df.columns.str.lstrip()

df.info()

for col in df.columns:
	if col in ('step_number','timestamp','pstate'):
		continue

	print(str(col) + ":")
	print("  min: " + str(df[col].min()))
	print("  avg: " + str(df[col].mean()))
	print("  max: " + str(df[col].max()))
	print("  med: " + str(df[col].median()))
	print("  std: " + str(df[col].std()))

plt.style.use(['dark_background'])

mpl.rcParams['axes.edgecolor'] = (0.600, 0.600, 0.600)

mpl.rcParams['axes.labelcolor'] = (0.800, 0.800, 0.800)
mpl.rcParams['axes.titlecolor'] = (0.800, 0.800, 0.800)

mpl.rcParams['text.color'] = (0.800, 0.800, 0.800)

mpl.rcParams['xtick.color'] = (0.600, 0.600, 0.600)
mpl.rcParams['ytick.color'] = (0.600, 0.600, 0.600)

mpl.rcParams['xtick.labelcolor'] = (0.800, 0.800, 0.800)
mpl.rcParams['ytick.labelcolor'] = (0.800, 0.800, 0.800)

if ('log_nvidia-smi' == log_name):
	df.plot(x='timestamp', y=['power.draw [W]'], figsize=(8,6))
	df.plot(x='timestamp', y=['fan.speed [%]'], figsize=(8,6))
else:
	df.plot(x='step_number', figsize=(8,6))

#fig, ax = plt.subplots()

#ax.plot(df['timestamp'], df['timestamp'])
#ax.plot(df['timestamp'], df['fan.speed [%]'])
#ax.plot(df['timestamp'], df['pstate'])
#ax.plot(df['timestamp'], df['utilization.gpu [%]'])
#ax.plot(df['timestamp'], df['utilization.memory [%]'])
#ax.plot(df['timestamp'], df['temperature.gpu'])
#ax.plot(df['timestamp'], df['power.draw [W]'])
#ax.plot(df['timestamp'], df['clocks.current.graphics [MHz]'])
#ax.plot(df['timestamp'], df['clocks.current.memory [MHz]'])

#ax.plot(df['timestamp'], df['potential_energy'])

#ax.plot(df['step_number'], df["drift"])
#ax.plot(df['step_number'], df["solver_run()"])
#ax.plot(df['step_number'], df["kick_1"])

plt.title(log_name)

plt.savefig(log_dir + '/' + log_name + '.png', dpi=300)

plt.show()
