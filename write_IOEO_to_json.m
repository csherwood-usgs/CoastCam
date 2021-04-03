% write_IOEO_to_dict - Convert Brittany's .mat format to json
clear
% specifiy input file
cal_fn = 'CACO02_C2_IOEOInitial_noT5.mat'

% make output file names
c = split(cal_fn, '.')
cal_fun_IO_json = [replace(c{1},'IOEO','IO'),'.json']
cal_fun_EO_json = [replace(c{1},'IOEO','EO'),'.json']
cal_fun_IOEO_json = [c{1},'.json']

% read calib file
load(cal_fn)

% write extrinsics in json format
otext = sprintf('{\n"x": %.3f,\n"y": %.3f,\n"z": %.3f,\n"a": %.3f,\n"t": %.3f,\n"r": %.3f\n}',...
    extrinsics)
fid = fopen(cal_fun_EO_json,'w');
fprintf(fid,'%s',otext);
fclose(fid);

% write intrinsics in json format
otext2 = sprintf('{\n"NU": %d,\n"NV": %d,\n"c0U": %f,\n"c0V": %f,\n"fx": %f,\n"fy": %f,\n"d1": %f,\n"d2": %f,\n"d3": %f,\n"t1": %f,\n"t2": %f\n}',...
    intrinsics)
fid = fopen(cal_fun_IO_json,'w');
fprintf(fid,'%s',otext2);
fclose(fid);

% combine both in json format
otext3 = sprintf('{\n"extrinsics": %s,\n"intrinsics": \n%s\n}',otext,otext2)
fid = fopen(cal_fun_IOEO_json,'w');
fprintf(fid,'%s',otext3);
fclose(fid);
