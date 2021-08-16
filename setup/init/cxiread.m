function [ cxi_struct ] = cxiread(filename)
%cxiread reads a CXI file resulting from a ptychographic reconstruction
%   Reads the entire HDF5 tree structure of a .cxi file
fid = H5F.open(filename,'H5F_ACC_RDONLY','H5P_DEFAULT');
cxi_struct = struct;
[status cxi_struct] = H5L.visit(fid,'H5_INDEX_NAME','H5_ITER_NATIVE',@iterFunc,cxi_struct);
end


function [status opdata_out] = iterFunc(group_id,name,opdata_in)
% group name
l_info = H5L.get_info(group_id,name,'H5P_DEFAULT');
if(l_info.type == H5ML.get_constant_value('H5L_TYPE_SOFT'))
    name = char(H5L.get_val(group_id,name,'H5P_DEFAULT'));
end
% get type and name
if(H5L.exists(group_id,name,'H5P_DEFAULT'))
    oid = H5O.open(group_id,name,'H5P_DEFAULT');    
    type = H5I.get_type(oid);   
    H5I.get_name(oid);
else
    type = -1;
end
if(type == H5ML.get_constant_value('H5I_DATASET'))
    data = H5D.read(oid,'H5ML_DEFAULT','H5S_ALL','H5S_ALL','H5P_DEFAULT');
    if(isstruct(data))
       f = fieldnames(data);
       if(numel(f) == 2 && f{1} == 'r' && f{2} == 'i')
           data = data.r + data.i*1i;
       end
    end
    opdata_out = SetProp(opdata_in,name,data);
    
else
    opdata_out = opdata_in;
end
status = 0;
end

function obj=SetProp(obj, prop, std_val, varargin)

prop = regexprep(prop, ' ','_');
S=struct('type','.','subs',regexp(prop,'\w*','match'));
obj = subsasgn(obj, S, (std_val));
end