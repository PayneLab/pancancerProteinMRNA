import cptac.pancan as pc
import sys, os

print(sys.argv[1])
currentdir = os.path.dirname(os.path.realpath('Make_Cancer_Delta_Corr_and_P_Value_Dataframe'))
parentdir = os.path.dirname(currentdir)
parentdir = os.path.dirname(parentdir)
sys.path.append(parentdir)


token = sys.argv[1] 


pc.download("pancanccrcc", box_token=token)
pc.download("pancanluad", box_token=token)
pc.download("pancanhnscc", box_token=token)
pc.download("pancanlscc", box_token=token)
pc.download("pancanucec", box_token=token)
