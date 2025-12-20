import pandas as pd
import matplotlib.pyplot as plt

plt.style.use(['dark_background'])

df = pd.read_csv('./out/timings_sim_step.csv')

df.plot(x="step_number")

#fig, ax = plt.subplots()

#ax.plot(df["step_number"], df["kick_0"])
#ax.plot(df["step_number"], df["drift"])
#ax.plot(df["step_number"], df["solver_run()"])
#ax.plot(df["step_number"], df["kick_1"])

plt.show()
