# NanoArray_Hopping
Hopping transport in NP array.
First
make sure all the .m files are in the same folder.
then
Start with RandRunLG_CM

This function defines some variables and calls the functions

RandlGMPGCM
and 
RandlLoopCM
(if a file has CM on the end then i modified it from YH's initial code)

it also has some code at the end to automatically save the output as figures and text.

first, just try to run it.

If it doesn't complete in 2 minutes, 
try smaller number for grids.
grids = 5
Trials = 1
NPgridrows = 6
you can also increase Vxmin to a larger voltage to help it run faster

You can also comment out all the sections at the bottom after the double for loops.

%%
output= 

you'll have to change the path variable if you do want to use the output saving.
