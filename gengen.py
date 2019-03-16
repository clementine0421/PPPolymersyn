import PDB_Synthesizer
import numpy as np
import randomgen


repeat = 1
entries = ["MB1","MB1","MB1","MB1","MB1",]*repeat
output_name= "BBBBB1"

#These settings work for intramonomer concatination (between MAs and MBs)
#link_setting = PDB_Synthesizer.LinkSetting(
#       connect_from="C",
#       connect_to="N",
#       #bond_offset=np.array([0., 0.5, -1.5])
#       #The above bond offset for MA0-MA0, MA1-MA1, MA0-MA1, MA1-MA0 concatination
#       bond_offset=np.array([-0.35, 0.6, 1.8])
#       #The above bond offset for MB0-MB0, MB1-MB1, MB0-MB1, MB1-MB0 concatination
#   )

link_setting_strings = ["{}{}".format(entries[i], entries[i+1]) for i in range(len(entries) - 1)]
link_settings = [randomgen.link_setting_map[lss] for lss in link_setting_strings]
print(link_setting_strings)

final_name = "./{}x{}.pdb".format(output_name,repeat)
PDB_Synthesizer.construct(
    entries,
    final_name,
    link_settings
)

print(final_name)
