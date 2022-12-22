import cptac.pancan as pc
import sys, os

print(sys.argv[1])
token = sys.argv[1] 


pc.download("pancanccrcc", box_token=token)
pc.download("pancanluad", box_token=token)
pc.download("pancanhnscc", box_token=token)
pc.download("pancanlscc", box_token=token)
pc.download("pancanucec", box_token=token)
