# %%
import numpy as np
import scipy.io as sio
import cv2 
import json
import matplotlib.pyplot as plt
import os
import sys

# %%
# %% parameters setting
subject = "Julie"
date = "20231012"
state = "IJVSmall"
with open("blood_vessel_segmentation_line_new.json") as f:
    file = json.load(f)
    length10mmEdge = file[subject][date]["length10mmEdge"]
    bound = file[subject][date]['bound']
    skin_x = file[subject][date]['skin']['x']
    fat_x = file[subject][date]['fat']['x']
    cca_v = file[subject][date][state]['cca']['v']

voxelLength = 0.25  # [mm]
gridNumIn10mm = int(10/voxelLength)

# %% [markdown]
# ## trim to tissue region ( delete not necessary region )

# %%
# %% main
# plot original image
image = cv2.imread(f'{subject}_{date}_{state}.png', cv2.IMREAD_GRAYSCALE)

## find tissue upper bound and lower bound
plt.imshow(image, cmap="gray")
plt.axhline(y=bound[0], color='red') # upperbound
plt.axhline(y=bound[1], color='red') # lowerbound
plt.plot()
plt.colorbar()
plt.title("original")
plt.show()
# plot tissue area only
plt.imshow(image[bound[0]:bound[1]], cmap="gray")
plt.colorbar()
plt.title("tissue (trim to skin surface)")
plt.show()
image = image[bound[0]:bound[1]]

# %% [markdown]
# ## Find skin and fat line

# %%
## find tissue skin,fat, thickness
# trim image to start from skin surface
plt.imshow(image, cmap="gray")
plt.axhline(y=skin_x, color='red') 
plt.axhline(y=fat_x, color='red') 
plt.plot([0, 0], length10mmEdge, "r.") # scale bar (10mm = ? pixels)
plt.colorbar()
plt.title("tissue (trim to skin surface)")
plt.show()

# %% [markdown]
# ## Using GUI to find IJV region

# %%
# plot the IJV bundary
img = image
im = []
XY = []
imgs = f'image_{subject}_{date}_{state}'
def left(ii,jj,l):
    for i in l:
        if (i[0]==ii and i[1]<jj):
            return True
    return False

def right(ii,jj,l):
    for i in l:
        if (i[0]==ii and i[1]>jj):
            return True
    return False

def up(ii,jj,l):
    for i in l:
        if (i[1]==jj and i[0]<ii):
            return True
    return False

def down(ii,jj,l):
    for i in l:
        if (i[1]==jj and i[0]>ii):
            return True
    return False

def drawing(event, x, y, flags, param):
    global img

    if event == 0 and flags == 1: #鼠標移動 and 左鍵按下
        if ([x,y] not in XY):
            XY.append([x,y])
        if ([x+1, y+1] not in XY):
            XY.append([x+1, y+1])
        if ([x-1, y-1] not in XY):
            XY.append([x-1, y-1])
        if ([x+1, y-1] not in XY):
            XY.append([x+1, y-1])
        if ([x-1, y+1] not in XY):
            XY.append([x-1, y+1])
        if ([x+1, y] not in XY):
            XY.append([x+1, y])
        if ([x-1, y] not in XY):
            XY.append([x-1, y])
        if ([x, y+1] not in XY):
            XY.append([x, y+1])
        if ([x, y-1] not in XY):
            XY.append([x, y-1])    
        cv2.circle(img, (x,y), 1, (71, 0, 255), -1)
        cv2.imshow('img1', img)

cv2.namedWindow("img1")
cv2.moveWindow("img1", 100, 100)
cv2.setMouseCallback("img1", drawing)
cv2.imshow("img1", img)

cv2.waitKey(0)
cv2.destroyAllWindows()


# # %%
# find the points(pixels) which are inside the bundary we plot
img_size = [image.shape[0], image.shape[1]]
if (XY):
    img2 = np.zeros(img_size, np.uint16)
    for m in XY:
        img2[m[1]][m[0]] = 255
    iii=0
    jjj=0
    start=0
    for i in img2:
        jjj=0
        if (right(iii,jjj,XY) or left(iii,jjj,XY) or up(iii,jjj,XY) or down(iii,jjj,XY)):
            start=1
        if (start):
            for j in i:
                if (right(iii,jjj,XY) and left(iii,jjj,XY) and up(iii,jjj,XY) and down(iii,jjj,XY)):
                    img2[jjj][iii]= 255
                jjj+=1
        start =0
        iii+=1
    plt.figure(figsize=(20,16))
    ax = plt.subplot(3,1,1)
    ax.imshow(img2)
    ax.axis("off")

# %%
# save ijv position
ijv_pos = np.where(img2==255)
plt.imshow(image, cmap="gray")
plt.plot(ijv_pos[1], ijv_pos[0], "r.", markersize=1)
plt.show()
np.save(f'{subject}_{date}_{state}_ijv_pos.npy', ijv_pos)

# %% [markdown]
# ## Find CCA region (assume CCA is a circle --> x^2 + y^2 <= r^2)

# %%
tissue = 'cca'
# highlight and catch target tissue
plt.imshow(image, cmap="gray")
plt.colorbar()
angles = np.linspace(0, 2*np.pi, num=100)
x = cca_v[2] * np.cos(angles) + cca_v[0]
y = cca_v[2] * np.sin(angles) + cca_v[1]
plt.plot(x, y, c="r")

# mesh coordinate
coordinates = np.meshgrid(
    np.arange(image.shape[0]), np.flip(np.arange(image.shape[1])))

# add dummy coordinate
coordinates = np.insert(coordinates, 0, 1, axis=0)

# sketch potential region
targetMatch = np.ones((image.shape[1], image.shape[0]), dtype=bool)
targetMatch = targetMatch & ((coordinates[1]-cca_v[1])**2 + (coordinates[2]-cca_v[0])**2 <= cca_v[2]**2)

# plot image (before scale)
legalRow, legalCol = np.where(targetMatch == True)
legalRow = targetMatch.shape[0]-legalRow
plt.imshow(image, cmap="gray")
plt.plot(legalRow, legalCol, "r.", markersize=3)
plt.title("catch {} (before scale)".format(tissue))
plt.show()

# save cca position
cca_pos = np.concatenate((legalRow.reshape(1,-1), legalCol.reshape(1,-1)), axis=0)
np.save(f'{subject}_{date}_{state}_cca_pos.npy', cca_pos)


