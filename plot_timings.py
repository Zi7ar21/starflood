import pandas as pd
import matplotlib.pyplot as plt

plt.style.use(['dark_background'])

df = pd.read_csv('./out/timings_sim_step.csv')
#df = pd.read_csv('./out/timings_vis_draw.csv')

for col in df.columns:
	if "step_number" == col:
		continue
	print(str(col) + ": ")
	print("  min: " + str(df[col].min()))
	print("  med: " + str(df[col].median()))
	print("  max: " + str(df[col].max()))
	print("  avg: " + str(df[col].mean()))
	print("  std: " + str(df[col].std()))

df.plot(x="step_number")

#fig, ax = plt.subplots()

#ax.plot(df["step_number"], df["kick_0"])
#ax.plot(df["step_number"], df["drift"])
#ax.plot(df["step_number"], df["solver_run()"])
#ax.plot(df["step_number"], df["kick_1"])

plt.show()
