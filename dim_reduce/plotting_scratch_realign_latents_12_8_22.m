% Cells files used can be downloaded from : https://drive.google.com/drive/folders/1-J6MjPZ7zlDePnbhOCeGRASkA8nIk0Ad?usp=sharing

file_path_18 = "C:\Users\abondy\Downloads\fulllatentdata_18.npz";
file_path_19 = "C:\Users\abondy\Downloads\fulllatentdata_19.npz";

latents_18 = readNPY(file_path_18);
latents_19 = readNPY(file_path_19);

Cells_18=load('C:\Users\abondy\Downloads\A297_2021_07_18_01_Cells.mat'); % you need thes for the violated field
Cells_19=load('C:\Users\abondy\Downloads\A297_2021_07_19_01_Cells.mat');

Cells_18 = add_first_click_state(Cells_18);
Cells_19 = add_first_click_state(Cells_19);

gpfa_18=format_latents(Cells_18,latents_18,'cpoke_in',[0.002 3.248],'trial_idx',~Cells_18.Trials.violated);
gpfa_19=format_latents(Cells_19,latents_19,'cpoke_in',[0.002 3.248],'trial_idx',~Cells_19.Trials.violated);

gpfa_18_first_click = realign_latents(gpfa_18,'first_click',[-0.2 1]);
gpfa_19_first_click = realign_latents(gpfa_19,'first_click',[-0.2 1]);

gpfa_18_cpoke_out = realign_latents(gpfa_18,'cpoke_out',[-1 1]);
gpfa_19_cpoke_out = realign_latents(gpfa_19,'cpoke_out',[-1 1]);

% permute back to trials x latents x time and save as npy
writeNPY(permute(gpfa_18_first_click.score,[2 3 1]),"C:\Users\abondy\Downloads\fulllatentdata_18_first_click.npy")
writeNPY(permute(gpfa_19_first_click.score,[2 3 1]),"C:\Users\abondy\Downloads\fulllatentdata_19_first_click.npy")

writeNPY(permute(gpfa_18_cpoke_out.score,[2 3 1]),"C:\Users\abondy\Downloads\fulllatentdata_18_cpoke_out.npy")
writeNPY(permute(gpfa_19_cpoke_out.score,[2 3 1]),"C:\Users\abondy\Downloads\fulllatentdata_19_cpoke_out.npy")