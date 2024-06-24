import pandas as pd
import os 

mol_close_contact = pd.read_csv('data_29900000.csv', names=range(1,201))
mol_close_contact.index +=1


faceon_pair = []
for mol_line in range(1,201):
    for mol_cul in range(1, 201):
        if(mol_close_contact.at[mol_line, mol_cul] >= 8):
            faceon_pair.append([mol_line, mol_cul])



#faceon_pair = {19:[5,83],71:[21,194]}

#for i,value in faceon_pair.items():
for i, value in enumerate(faceon_pair,start=1):
    os.system('echo RUNNING '+str(i))
    os.system('echo pair '+ str(value[0]) +' and'+ str(value[1]))
    os.system('echo running lumo')
    os.system('/home/21141213516/YWQ/calc_J -p_1 /home/21141213516/YWQ//n32990v//'+str(value[0])+'.pun -orb_ty_1 LUMO -p_2 /home/21141213516/YWQ//n32990v//'+str(value[1])+'.pun -orb_ty_2 LUMO -p_P /home/21141213516/YWQ//n32990v//'+str(i)+'_pair.pun')
    os.system('echo running homo')
    os.system('/home/21141213516/YWQ/calc_J -p_1 /home/21141213516/YWQ//n32990v//'+str(value[0])+'.pun -p_2 /home/21141213516/YWQ//n32990v//'+str(value[1])+'.pun -p_P /home/21141213516/YWQ//n32990v//'+str(i)+'_pair.pun')
    os.system('echo '+str(i)+' has finished')
