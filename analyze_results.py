import pandas as pd
import matplotlib.pyplot as plt

#filename = "./log_nvidia-smi_2025-12-20.csv"
filename = './out/statistics.csv'
#filename = './out/timings_sim_step.csv'
#filename = './out/timings_vis_draw.csv'

df = pd.read_csv(filename)

for col in df.columns:
	if "step_number" == col:
		continue
	print(str(col) + ": ")
	print("  min: " + str(df[col].min()))
	print("  med: " + str(df[col].median()))
	print("  max: " + str(df[col].max()))
	print("  avg: " + str(df[col].mean()))
	print("  std: " + str(df[col].std()))

plt.style.use(['dark_background'])

df.plot(x="step_number")

#fig, ax = plt.subplots()

#ax.plot(df["step_number"], df["kick_0"])
#ax.plot(df["step_number"], df["drift"])
#ax.plot(df["step_number"], df["solver_run()"])
#ax.plot(df["step_number"], df["kick_1"])

plt.title(filename)

plt.savefig('./statistics.png', dpi=300)

plt.show()
