# ABANICCO

![ABANICCO_LOGO](https://user-images.githubusercontent.com/49556605/201992997-359d29d3-8c88-4d22-8bba-1517db572d01.png)

This is the first version of ABANICCO: an AB ANgular Illustrative Classification of COlor. The corresponding paper can be found in [SENSORS](https://www.mdpi.com/1424-8220/23/6/3338)

This version requires MATLAB2021 or higher. You can also [![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=lauranicolass/ABANICCO)

To use it, just run examplerun.mat. All instructions, examples and options are present in that file.


Different applications will need different actions:

## For **MULTIPLEX** images:

We recommend running ABANICCO on the separate channels and then comparing that with the merged image. The easiest way of doing so is a tiles image of the channels. However, the best way of doing it is individually. 

There is no need to pre-process the images. Increasing brightness and contrast will make visualization clearer, but will also amplify noise:

*Multiplex channels:*
![m_dark](https://user-images.githubusercontent.com/49556605/204319951-b3b22901-e5cc-4432-88d6-756a2848163d.png)

*Multiplex channels with increased brightness and contrast:*
![m_bright](https://user-images.githubusercontent.com/49556605/204319941-00d2eec0-fdf3-4f33-89b0-1dd6f8ece99c.png)

