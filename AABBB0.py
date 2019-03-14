import PDB_Synthesizer
import numpy as np

entries = ["MB1"]*6
repeat = 6


#These settings work for  MA0-MA0
link_setting = PDB_Synthesizer.LinkSetting(
       connect_from="C",
       connect_to="N",
       #bond_offset=np.array([0., 0.5, -1.5])
       #The above bond offset for MA0, MA1 concatination
       bond_offset=np.array([-0.1, 1, 1.6])
       #The above bond offset for MB0, MB1 concatination
   )



PDB_Synthesizer.construct(
    entries,
    "./{0}x{1}.pdb".format("MB1", repeat),
    [link_setting]*(repeat-1)
)
