# Holistic Engine Logic Program
Welcome to HELP! (Holistic Engine Logic Program), where "help" represents the anguished final cry of a struggling engineer before he or she disappears beneath the waves of an unholy ammount of STEM coursework. 

This is an open source MATLAB program that interfaces with NASA CEA to size and analyze most aspects of liquid rocket engine design. Specific features include (some to a limited extent or still WIP): engine geometry, contour export to CAD, injector sizing, regenerative cooling channel sizing, thermal analysis, structural analysis, fluid feed system sizing.

![image](https://user-images.githubusercontent.com/24968256/199154329-5501d4da-2fab-47af-bade-da9d5188470c.png)

![750 lbf regen](https://user-images.githubusercontent.com/24968256/201764015-f81a08c4-ce96-4222-85fe-85c37f19cd4a.png)

![750 lbf regen 2](https://user-images.githubusercontent.com/24968256/201764028-78e5458a-d8ce-4005-8347-28a7fe914668.PNG)

## Operation
1. Using the imperial or metric examples in /input, create a input excel file with your given run parameters.
2. Run main.m in the main directory after changing line 27 to your input file name.
3. Your program outputs will be outputed into the /output folder. (currently just outputs to console)

## Dependencies
- The REFPROP matlab wrapper is not required for general engine sizing but is required for fluid system, regenerative cooling, and injector sizing. Due to copyright rules, I cannot include REFPROP as part of this repository. enableREFPROP can be set to false in your input.xslx to prevent errors when you do not have it downloaded.

## Support
If you need any further help or would like to be added as a collaborator, email kamon.blong@gmail.com

