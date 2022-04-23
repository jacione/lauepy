import json
import numpy as np
print('################### GRAIN DICT #######################')
with open('/home/34idc-work/2022/LauePUP422/Analysis/lauepy_output/scan_0014/grains/grains.json') as f:
    grains = json.load(f)

print(f'{"Grain":>12}{"RMS":>10}{"Peaks":>10}{"Dist":>10}{"Pstdev":>10}{"Frames":>10}')
for grain in grains:
    g = grains[grain]
    pos_std = np.mean(np.std(g['Positions'], axis=0))
    num_frames = len(g['Frames'])
    if num_frames > 1 and pos_std < 1:
        goodness = np.log(num_frames * np.sqrt(g['Avg_Peaks']) / g['Avg_RMS'] / pos_std)
        if goodness > 5:
            print(f"{grain:>12}"
                  f"{g['Avg_RMS']:>10}"
                  f"{g['Avg_Peaks']:>10}"
                  f"{round(g['Avg_Dist'], 1):>10}"
                  f"{np.around(pos_std, 2):>10}"
                  f"{num_frames:>10}"
                  f"{np.around(goodness, 5):>11}"
                  )
