import cv2
import matplotlib.pyplot as plt
import pandas as pd
import json
cv2.startWindowThread()


img = cv2.imread('TileScan 3_Region3_HFD_No_AOMDSS_Merged_RAW_ch00.tif')
# img = cv2.imread('V1_Adult_Mouse_Brain_image.tif')
columns = ["barcode","in_tissue","array_row","array_col","pxl_row_in_fullres", "pxl_col_in_fullres"]
df_locs = pd.read_csv("tissue_positions_list.csv", header=None)
df_locs.columns = columns
json_file = open("scalefactors_json.json", "r")
scale_dict = data = json.loads(json_file.read())
radius = round(scale_dict["spot_diameter_fullres"]/2)+1


for ind, row in df_locs.iterrows():
    # print(row)
    if row["in_tissue"]==1:
    
        x_coor = row["pxl_col_in_fullres"]
        y_coor = row["pxl_row_in_fullres"]
        # cv2.circle(img,(x_coor,y_coor), radius, (0,0,255), -1)
        barcode= row["barcode"]
        cropped = img[y_coor-1-(radius):y_coor+1+(radius),  x_coor-1-(radius):x_coor+(radius)+1]
        cv2.imwrite(f"./tiles/{barcode}.tif", cropped)

		



# cv2.imshow('imag2e',img)
cv2.imwrite("deneme.tif", img)
cv2.waitKey(1000)