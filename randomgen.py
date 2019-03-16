"""
    random.py allows construcing PP polymer with random assignment of MAs and MBs
"""
import PDB_Synthesizer
import numpy as np
import random

entries = []
output_name = "Random16"

for i in range(16):
    entries.append(random.choice(["MA0", "MB1","MB0","MA1"]))



link_setting_map = {
    "MA0MA0": PDB_Synthesizer.LinkSetting(connect_from="C", connect_to="N", bond_offset=np.array([0., 0.5, -1.5])),
    "MA1MA0": PDB_Synthesizer.LinkSetting(connect_from="C", connect_to="N", bond_offset=np.array([0., 0.5, -1.5])),
    "MA1MA1": PDB_Synthesizer.LinkSetting(connect_from="C", connect_to="N", bond_offset=np.array([0., 0.5, -1.5])),
    "MB1MB1": PDB_Synthesizer.LinkSetting(connect_from="C", connect_to="N", bond_offset=np.array([-0.3,-1, 1.6])),
    "MB0MB0": PDB_Synthesizer.LinkSetting(connect_from="C", connect_to="N", bond_offset=np.array([-0.3,-1, 1.6])),
    "MB0MB1": PDB_Synthesizer.LinkSetting(connect_from="C", connect_to="N", bond_offset=np.array([-0.3,-1, 1.6])),
    "MA0MB0": PDB_Synthesizer.LinkSetting(connect_from="C", connect_to="N", bond_offset=np.array([-0.2,-0.8,-1.7])),
    "MA0MB1": PDB_Synthesizer.LinkSetting(connect_from="C", connect_to="N", bond_offset=np.array([-0.2,-0.8,-1.7])),
    "MA1MB1": PDB_Synthesizer.LinkSetting(connect_from="C", connect_to="N", bond_offset=np.array([-0.2,-0.8,-1.7])),
    "MA1MB0": PDB_Synthesizer.LinkSetting(connect_from="C", connect_to="N", bond_offset=np.array([-0.2,-0.8,-1.7])),
    "MB1MB0": PDB_Synthesizer.LinkSetting(connect_from="C", connect_to="N", bond_offset=np.array([-0.3,-1, 1.6])),
    "MB1MB0": PDB_Synthesizer.LinkSetting(connect_from="C", connect_to="N", bond_offset=np.array([-0.3,-1, 1.6])),
    "MB0MA1": PDB_Synthesizer.LinkSetting(connect_from="C", connect_to="N", bond_offset=np.array([0.2, -0.6, 1.6])),
    "MB0MA0": PDB_Synthesizer.LinkSetting(connect_from="C", connect_to="N", bond_offset=np.array([0.2, -0.6, 1.6])),
    "MB1MA0": PDB_Synthesizer.LinkSetting(connect_from="C", connect_to="N", bond_offset=np.array([0.2, -0.6, 1.6])),
    "MB1MA1": PDB_Synthesizer.LinkSetting(connect_from="C", connect_to="N", bond_offset=np.array([0.2, -0.6, 1.6]))
}

if __name__ == "__main__":

    link_setting_strings = ["{}{}".format(entries[i], entries[i+1]) for i in range(len(entries) - 1)]
    link_settings = [link_setting_map[lss] for lss in link_setting_strings]
    final_name = "./{}.pdb".format(output_name)
    print(final_name)


    PDB_Synthesizer.construct(
        entries,
        #"./{0}x{1}.pdb".format("MB1", repeat),
        "./{}.pdb".format(output_name),
        link_settings
        #[link_setting]*(repeat-1)
    )
