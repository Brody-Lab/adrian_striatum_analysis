function log = generate_implant_log()

    i=0;
    P=get_parameters;
    
    % look up table defining each probe implantation.
    % to add a new entry just copy paste one of the existing entries to
    % the bottom and replace fields as appropriate.
    

    %% ML angle is positive when tip is more lateral than insertion site
    %% AP angle is positive when tip is more anterior than insertion site
    %% shank_orientation = 0 when parallel to sagittal plane, 90 when parallel to coronal plane. values are accurate to 5 degrees unless more precise measurements were taken (as stated when applicable)

    %% A242
    i=i+1;
    log(i).rat = "A242";
    log(i).date = "2019-05-20";    
    log(i).serial = "17131311621";
    log(i).probe_type="10";
    log(i).insertion_depth_mm = 7.6;    
    log(i).ML_mm = 4.95;
    log(i).AP_mm = -2.35;
    log(i).coronal_angle_deg = 5; % positive means lateral as you go down
    log(i).sagittal_angle_deg = 0; % positive means you go anterior as you go down
    log(i).hemisphere = "right";
    log(i).shank_orientation = 0;    
    log(i).ap_group = 4;    
    log(i).region(1).name = "Amyg";
    log(i).region(1).electrode_range = [1 70];
    log(i).region(2).name = "TS";
    log(i).region(2).electrode_range = [71 384];
    log(i).region(3).name = "S1";
    log(i).region(3).electrode_range = [385 960];
    log(i).comments="";
    

    %% A230
    i=i+1;
    log(i).serial = "18005106831";
    log(i).probe_type="10";    
    log(i).ap_group = 1;
    log(i).date = "2019-07-02";
    log(i).ML_mm = 0.5;
    log(i).AP_mm = 4.0;
    log(i).insertion_depth_mm = 7.5;
    log(i).coronal_angle_deg = 26;
    log(i).sagittal_angle_deg = -29;
    log(i).hemisphere = "right";
    log(i).rat = "A230";
    log(i).shank_orientation = 45;
    log(i).comments = "haven't gone through the histology yet but it is available";

    i=i+1;
    log(i).serial = "17131308571";
    log(i).probe_type="10";    
    log(i).ap_group = 3;
    log(i).date = "2019-07-02";
    log(i).ML_mm = 4.0;
    log(i).AP_mm = -0.8;
    log(i).insertion_depth_mm = 6.6;
    log(i).coronal_angle_deg = -2; % positive means lateral as you go down
    log(i).sagittal_angle_deg = 0; % positive means you go anterior as you go down
    log(i).hemisphere = "left";
    log(i).rat = "A230";
    log(i).region(1).name = "IC"; 
    log(i).region(1).electrode_range = [1 200];
    log(i).region(2).name = "Str3";
    log(i).region(2).electrode_range = [201 470];
    log(i).region(3).name = "S1";
    log(i).region(3).electrode_range = [471 960];
    log(i).comments = "confirmed in ephys that waveforms in region 1 look a lot like axons  :(. I really don't think it's GP.";

    i=i+1;
    log(i).serial = "17131308411";
    log(i).probe_type="10";    
    log(i).ap_group=4;
    log(i).date = "2019-07-02";
    log(i).ML_mm = 5.0;
    log(i).AP_mm = -2.2;
    log(i).insertion_depth_mm = 8.6;
    log(i).coronal_angle_deg = 5;
    log(i).sagittal_angle_deg = 0;
    log(i).hemisphere = "right";
    log(i).rat = "A230";
    log(i).region(1).name = "PCx"; % piriform cortex
    log(i).region(1).electrode_range = [1 80];
    log(i).region(2).name = "DEN"; % dorsal endopeduncular nucleus
    log(i).region(2).electrode_range = [81 170];
    log(i).region(3).name = "Amyg";
    log(i).region(3).electrode_range = [171 230];
    log(i).region(4).name = "TS";
    log(i).region(4).electrode_range = [231 480];
    log(i).region(5).name = "S1";
    log(i).region(5).electrode_range = [481 960];
    log(i).shank_orientation = 0;
    log(i).comments="I think I've underestimated the distances from the bottom on this probe. Could have TS cells at slightly higher electrodes. I think this because the probe was inserted longer than I account for brain areas looking at the histology.";

    %% A241
    i=i+1;
    log(i).serial = "18194823631";
    log(i).probe_type="10";    
    log(i).ap_group=2;
    log(i).date = "2019-09-11";
    log(i).ML_mm = 2.15;
    log(i).AP_mm = 0.7; 
    log(i).insertion_depth_mm = 10; 
    log(i).coronal_angle_deg = 2; % positive means lateral as you go down
    log(i).sagittal_angle_deg = 0; % positive means you go anterior as you go down
    log(i).hemisphere = "right";
    log(i).rat = "A241";
    log(i).region(1).name = "PO"; 
    log(i).region(1).electrode_range = [1 300];
    log(i).region(2).name = "BNST"; 
    log(i).region(2).electrode_range = [301 520];
    log(i).region(3).name = "Str2";
    log(i).region(3).electrode_range = [521 740];
    log(i).region(4).name = "M1";
    log(i).region(4).electrode_range = [741 960];
    log(i).shank_orientation = 0;
    log(i).comments = "from histology, AP looks more like -0.4. all electrodes in the brain, so depth must be >9.6. the VP, striatum border is extremely hard to estimate until i do histology";

    i=i+1;
    log(i).serial = "18194823302";
    log(i).probe_type="10";    
    log(i).ap_group=3;
    log(i).date = "2019-09-11";
    log(i).ML_mm = 4.0;
    log(i).AP_mm = -0.6; % from histology, looks more like -1.3
    log(i).insertion_depth_mm = 10; % all probes are in, so >9.8
    log(i).coronal_angle_deg = -2; % positive means lateral as you go down
    log(i).sagittal_angle_deg = 0; % positive means you go anterior as you go down
    log(i).hemisphere = "left";
    log(i).rat = "A241";
    log(i).region(1).name = "Amyg";
    log(i).region(1).electrode_range = [1 209];    
    log(i).region(2).name = "IC";
    log(i).region(2).electrode_range = [210 509];    
    log(i).region(3).name = "Str3";
    log(i).region(3).electrode_range = [510 710];
    log(i).region(4).name = "S1";
    log(i).region(4).electrode_range = [710 960];
    log(i).shank_orientation = 0;
    log(i).comments = "from histology, AP looks more like -1.3. all electrodes in the brain, so depth must be >9.6. ";
    

    %% A243
    i=i+1;
    log(i).serial = "18194824132";
    log(i).probe_type="10";    
    log(i).ap_group=2;
    log(i).date = "2019-09-13";
    log(i).ML_mm = 2.45;
    log(i).AP_mm = 0.7; % 
    log(i).insertion_depth_mm = 9;
    log(i).coronal_angle_deg = 0; 
    log(i).sagittal_angle_deg = 0; 
    log(i).hemisphere = "right";
    log(i).rat = "A243";
    log(i).region(1).name = "Amyg";
    log(i).region(1).electrode_range = [1 69];
    log(i).region(2).name = "PO";
    log(i).region(2).electrode_range = [70 169];    
    log(i).region(3).name = "VP";
    log(i).region(3).electrode_range = [170 269];
    log(i).region(4).name = "GP";
    log(i).region(4).electrode_range = [270 419];
    log(i).region(5).name = "Str2";
    log(i).region(5).electrode_range = [420 670];
    log(i).region(6).name = "M1";
    log(i).region(6).electrode_range = [720 920];
    log(i).shank_orientation = 0;
    log(i).comments="from histology, looks more like -0.5 AP.  depth based on fact that bank 2, channel 112 is most superficial probe in brain.";

    i=i+1;
    log(i).serial = "18194823211";
    log(i).probe_type="10";    
    log(i).ap_group=3;
    log(i).date = "2019-09-13";
    log(i).ML_mm = 4.0;
    log(i).AP_mm = -0.6; % 
    log(i).insertion_depth_mm = 9.45; % 
    log(i).coronal_angle_deg = -2; 
    log(i).sagittal_angle_deg = 0; 
    log(i).hemisphere = "left";
    log(i).rat = "A243";
    log(i).region(1).name = "Amyg";
    log(i).region(1).electrode_range = [1 149];    
    log(i).region(2).name = "GP";
    log(i).region(2).electrode_range = [150 449];
    log(i).region(3).name = "Str3";
    log(i).region(3).electrode_range = [450 674];
    log(i).region(4).name = "S1";
    log(i).region(4).electrode_range = [675 925];    
    log(i).shank_orientation = 0;
    log(i).comments="from histology, looks more like -1.2 AP. depth based on that fact bank 2, channel 157 is most superficial probe in brain";

    % X046
    i=i+1;
    log(i).rat = "X046";
    log(i).serial = "18194823122";
    log(i).probe_type="10";    
    log(i).ap_group=1;
    log(i).date = "2019-09-15";
    log(i).ML_mm = 1.3; % ~1.2 (still in M2)
    log(i).AP_mm = 1.9; % ~1.5 from histology
    log(i).insertion_depth_mm = 7.4;
    log(i).coronal_angle_deg = 10; 
    log(i).sagittal_angle_deg = 0; 
    log(i).hemisphere = "left";
    log(i).region(1).name = "NAc";
    log(i).region(1).electrode_range = [1 149];
    log(i).region(2).name = "ADS";
    log(i).region(2).electrode_range = [150 456];
    log(i).region(3).name = "FOF";
    log(i).region(3).electrode_range = [457 960];
    log(i).comments="";    
    log(i).shank_orientation = [];

    % A249 - no histology
    i=i+1;
    log(i).serial = "18194819132";
    log(i).probe_type="10";    
    log(i).ap_group=1;
    log(i).date = "2020-02-04";
    log(i).ML_mm = 2.1;
    log(i).AP_mm = 2.2;
    log(i).insertion_depth_mm = 6.8;
    log(i).coronal_angle_deg = 0; % positive means lateral as you go down
    log(i).sagittal_angle_deg = -5; % positive means you go anterior as you go down
    log(i).hemisphere = "right";
    log(i).rat = "A249";
    log(i).region(1).name = "NAc";
    log(i).region(1).electrode_range = [1 125];    
    log(i).region(2).name = "ADS";
    log(i).region(2).electrode_range = [126 420];
    log(i).region(3).name = "M2";
    log(i).region(3).electrode_range = [421 960];    
    log(i).shank_orientation = 90;
    log(i).comments = "haven't looked yet at histology";
    
    i=i+1;
    log(i).serial = "18194819321";
    log(i).probe_type="10";    
    log(i).ap_group=1;
    log(i).date = "2020-09-01";
    log(i).ML_mm = 3.1;
    log(i).AP_mm = 2.2;
    log(i).insertion_depth_mm = 6.8;
    log(i).coronal_angle_deg = 0; % positive means lateral as you go down
    log(i).sagittal_angle_deg = -5; % positive means you go anterior as you go down
    log(i).hemisphere = "right";
    log(i).rat = "A256";    
    log(i).region(1).name = "NAc";
    log(i).region(1).electrode_range = [1 100];    
    log(i).region(2).name = "ADS";
    log(i).region(2).electrode_range = [101 420];
    log(i).region(3).name = "M1";
    log(i).region(3).electrode_range = [421 960];    
    log(i).shank_orientation = 90;    
    log(i).comments = "haven't looked yet at histology";
    
    
    
    % A297 - TS
    i=i+1;
    log(i).serial = "19076606841";
    log(i).probe_type="10";    
    log(i).ap_group=4;
    log(i).date = "2021-07-14";
    log(i).ML_mm = 5.0;
    log(i).AP_mm = -2.1;
    log(i).insertion_depth_mm = 7.0;
    log(i).coronal_angle_deg = 5; % positive means lateral as you go down
    log(i).sagittal_angle_deg = 0; % positive means you go anterior as you go down
    log(i).hemisphere = "right";
    log(i).rat = "A297";    
    log(i).region(1).name = "Amyg";
    log(i).region(1).electrode_range = [1 60];    
    log(i).region(2).name = "TS";
    log(i).region(2).electrode_range = [61 420];
    log(i).region(3).name = "S1";
    log(i).region(3).electrode_range = [421 960];
    log(i).comments="";
    log(i).shank_orientation = 90;    
    
    % A297 - ADS
    i=i+1;
    log(i).serial = "18194819321";
    log(i).probe_type="10";    
    log(i).ap_group=1;
    log(i).date = "2021-07-14";
    log(i).ML_mm = 1.0;
    log(i).AP_mm = 2.0;
    log(i).insertion_depth_mm = 7;
    log(i).coronal_angle_deg = 12; % positive means lateral as you go down
    log(i).sagittal_angle_deg = 0; % positive means you go anterior as you go down
    log(i).hemisphere = "right";
    log(i).rat = "A297";    
    log(i).region(1).name = "ADS";
    log(i).region(1).electrode_range = [1 420];    
    log(i).region(2).name = "FOF";
    log(i).region(2).electrode_range = [421 960];
    log(i).comments="";
    log(i).shank_orientation = 90;        
    
    % A294 - TS
    i=i+1;
    log(i).serial = "19076606882";
    log(i).probe_type="10";    
    log(i).ap_group=4;
    log(i).date = "2021-06-09";
    log(i).ML_mm = 5.0;
    log(i).AP_mm = -2.1;
    log(i).insertion_depth_mm = 7.0;
    log(i).coronal_angle_deg = 5; % positive means lateral as you go down
    log(i).sagittal_angle_deg = 0; % positive means you go anterior as you go down
    log(i).hemisphere = "right";
    log(i).rat = "A294";    
    log(i).region(1).name = "Amyg";
    log(i).region(1).electrode_range = [1 60];    
    log(i).region(2).name = "TS";
    log(i).region(2).electrode_range = [61 420];
    log(i).region(3).name = "S1";
    log(i).region(3).electrode_range = [421 960];
    log(i).shank_orientation = 90;    
    log(i).comments="";    
    
    % A294 - ADS
    i=i+1;
    log(i).serial = "18194819132";
    log(i).probe_type="10";    
    log(i).ap_group=1;
    log(i).date = "2021-06-09";
    log(i).ML_mm = 1.0;
    log(i).AP_mm = 2.0;
    log(i).insertion_depth_mm = 7;
    log(i).coronal_angle_deg = 12; % positive means lateral as you go down
    log(i).sagittal_angle_deg = 0; % positive means you go anterior as you go down
    log(i).hemisphere = "right";
    log(i).rat = "A294";    
    log(i).region(1).name = "ADS";
    log(i).region(1).electrode_range = [1 420];    
    log(i).region(2).name = "FOF";
    log(i).region(2).electrode_range = [421 960];
    log(i).comments="";
    log(i).shank_orientation = 90;   
    
    % A295
    i=i+1;
    log(i).serial = "18194819132";
    log(i).probe_type="10";    
    log(i).ap_group=3;
    log(i).date = "2021-09-29";
    log(i).ML_mm = 4.0;
    log(i).AP_mm = -0.4;
    log(i).insertion_depth_mm = 7.0;
    log(i).coronal_angle_deg = -2; % positive means lateral as you go down
    log(i).sagittal_angle_deg = 0; % positive means you go anterior as you go down
    log(i).hemisphere = "right";
    log(i).rat = "A295";    
    log(i).region(1).name = "Str3";
    log(i).region(1).electrode_range = [200 425];    
    log(i).region(2).name = "M2";
    log(i).region(2).electrode_range = [426 960];
    log(i).shank_orientation = 90;    
    log(i).comments = "histology needs to be examined";
    
    % A295 
    i=i+1;
    log(i).serial = "19076606882";
    log(i).probe_type="10";    
    log(i).ap_group=3;
    log(i).date = "2021-09-29";
    log(i).ML_mm = 4.0;
    log(i).AP_mm = -0.4;
    log(i).insertion_depth_mm = 7;
    log(i).coronal_angle_deg = -2; % positive means lateral as you go down
    log(i).sagittal_angle_deg = 0; % positive means you go anterior as you go down
    log(i).hemisphere = "left";
    log(i).rat = "A295";    
    log(i).region(1).name = "Str3";
    log(i).region(1).electrode_range = [200 425];    
    log(i).region(2).name = "M2";
    log(i).region(2).electrode_range = [426 960];
    log(i).shank_orientation = 90;     
    log(i).comments = "histology needs to be examined";    
    
    % T219
    i=i+1;
    log(i).rat = "T219";
    log(i).serial = "18194824092";
    log(i).probe_type="10";    
    log(i).ap_group=1;
    log(i).date = "2019-11-17";
    log(i).ML_mm = 2.4;
    log(i).AP_mm = 1;
    log(i).insertion_depth_mm = 8;
    log(i).coronal_angle_deg = 0; % positive means lateral as you go down
    log(i).sagittal_angle_deg = 15; % positive means you go anterior as you go down
    log(i).hemisphere = "left";
    log(i).region(1).name = "NAc";
    log(i).region(1).electrode_range = [1 261];
    log(i).region(2).name = "ADS";
    log(i).region(2).electrode_range = [262 532];
    log(i).region(3).name = "M1";
    log(i).region(3).electrode_range = [532 960];
    log(i).comments="";
    log(i).shank_orientation = [];
    
     % A327, left ADS (based on A324)
    i=i+1;
    log(i).rat = "A327";
    log(i).serial = "19076606842";
    log(i).probe_type="10";    
    log(i).ap_group=1;
    log(i).date = "2023-08-08";
    log(i).ML_mm = 3;
    log(i).AP_mm = 1.9;
    log(i).insertion_depth_mm = 7.9;
    log(i).coronal_angle_deg = -10;
    log(i).sagittal_angle_deg = 0;
    log(i).hemisphere = "left";
    log(i).region(1).name = "NAc";
    log(i).region(1).electrode_range = [1 220];
    log(i).region(2).name = "ADS";
    log(i).region(2).electrode_range = [221 530];
    log(i).region(3).name = "M1";
    log(i).region(3).electrode_range = [610 960];
    log(i).shank_orientation = [];    
    log(i).comments = "regions based on A324 extraplotion prior to having histology";
    
    
    % A327, right ADS (based on A324)
    i=i+1;
    log(i).rat = "A327";
    log(i).serial = "19051016061";
    log(i).probe_type="10";    
    log(i).ap_group=1;
    log(i).date = "2023-08-08";
    log(i).ML_mm = 3;
    log(i).AP_mm = 1.9;
    log(i).insertion_depth_mm = 7.6;
    log(i).coronal_angle_deg = -10;
    log(i).sagittal_angle_deg = 0;
    log(i).hemisphere = "right";
    log(i).region(1).name = "NAc";
    log(i).region(1).electrode_range = [1 190];
    log(i).region(2).name = "ADS";
    log(i).region(2).electrode_range = [191 500];
    log(i).region(3).name = "M1";
    log(i).region(3).electrode_range = [580 960];
    log(i).shank_orientation = [];     
    log(i).comments = "regions based on A324 extrapolation prior to having histology"   
    
     % A327, left TS  (estimated from A325 histology and A327 ephys)
    i=i+1;
    log(i).rat = "A327";
    log(i).serial = "19076610451";
    log(i).probe_type="10";    
    log(i).ap_group=4;
    log(i).date = "2023-08-08";
    log(i).ML_mm = 5;
    log(i).AP_mm = -2;
    log(i).insertion_depth_mm = 7.56;
    log(i).coronal_angle_deg = 5;
    log(i).sagittal_angle_deg = 0;
    log(i).hemisphere = "left";
%     log(i).region(1).name = "Amyg";
%     log(i).region(1).electrode_range = [1 100];   
%     log(i).region(2).name = "TS";
%     log(i).region(2).electrode_range = [101 460];
%     log(i).region(3).name = "S1";
%     log(i).region(3).electrode_range = [461 960];   
%     log(i).shank_orientation = [];       
    
    
    log(i).region(1).name = "TS";
    log(i).region(1).electrode_range = [51 450];
    log(i).region(2).name = "S1";
    log(i).region(2).electrode_range = [491 960];
    log(i).shank_orientation = [];   
    
    log(i).comments = "regions based on A324 histology and A327 phys";        
 
    
    % A327, right TS (estimated from A325 histology and A327 ephys)
    i=i+1;
    log(i).rat = "A327";
    log(i).serial = "19076604731";
    log(i).probe_type="10";    
    log(i).ap_group=4;
    log(i).date = "2023-08-08";
    log(i).ML_mm = 5;
    log(i).AP_mm = -2.1;
    log(i).insertion_depth_mm = 7.7; 
    log(i).coronal_angle_deg = 5;
    log(i).sagittal_angle_deg = 0;
    log(i).hemisphere = "right";
    %log(i).region(1).name = "Amyg";
    %log(i).region(1).electrode_range = [1 110];   
    log(i).region(1).name = "TS";
    log(i).region(1).electrode_range = [50 430];
    log(i).region(2).name = "S1";
    log(i).region(2).electrode_range = [471 960];
    log(i).shank_orientation = [];     
    log(i).comments = "regions based on A324 histology and A327 phys";    
        
%%  write implant log as json file
    fid=fopen(P.implant_log_path,'w');
    fprintf(fid,jsonencode(log,'PrettyPrint',true));
    fclose(fid);
 
end
