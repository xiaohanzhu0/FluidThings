function [data, vartype, varname, ncvs, gridfile] = readTurtleFields(filename)
%
%  Read files containing Turtle fields, created by calls to Turtle::openAndWriteHeader()
%  and multiple calls to Turtle::writeSingleField().
%  Can also read files created by derived classes, by skipping parts of the header.
%
%clear; filename = 'snapshots/blk5_074_sgs_n0067237_t1.100e+02.res.avg_xz';

%  Parameters from Turtle_parameters.h
MAGICNUMBER    = 751123 ;
IO_NAME_LENGTH =    200 ;

%  Parameters from VariableList.h
CV_SCALAR =  1 ;
CV_VECTOR =  2 ;
CV_TENSOR =  3 ;
FA_SCALAR = 11 ;
FA_VECTOR = 12 ;
FA_TENSOR = 13 ;
NO_SCALAR = 21 ;
NO_VECTOR = 22 ;
NO_TENSOR = 23 ;


%  Open and check big-endian
endian = 'ieee-le';
fid = fopen(filename, 'r', endian);
magicnbr = fread(fid, 1, 'integer*4');
fclose(fid);
if ( magicnbr ~= MAGICNUMBER )
  endian = 'ieee-be';
end
fid = fopen(filename, 'r', endian);
magicnbr = fread(fid, 1, 'integer*4');
if ( magicnbr ~= MAGICNUMBER )
  disp('Error: can not determine the endian...')
  return;
end

%  Read the Turtle header
version = fread(fid, 1, 'integer*4');
if version < 10000*(2) + 100*(0) + (36),
  disp('Error: this version only reads Tortuga2 files')
  return
end
nbrBlocks = fread(fid, 1, 'integer*4');
gridfile = fread(fid, IO_NAME_LENGTH)';
i = find( gridfile == 0 ,1);
gridfile = char( gridfile(1:i-1) );
clear i
ncvs = reshape( fread(fid, 3*nbrBlocks, 'integer*4') , 3, nbrBlocks)';
remainingHeaderSize = fread(fid,1,'integer*4');

%  Skip remaining header (possibly created by derived class)
if remainingHeaderSize > 0,
  %disp(['  skipping ',int2str(remainingHeaderSize),' bytes in header...'])
  dummy = fread(fid, remainingHeaderSize, 'uchar');  clear dummy
end

%  Keep reading while there is still data
v=0;
tmp = fread(fid, 1, 'integer*4');   %  try reading the next vartype
while feof(fid)==0,
  v=v+1;
  vartype(v) = tmp;

  %  Read the rest of the variable header
  precision = fread(fid, 1, 'integer*4');
  dummy = fread(fid, IO_NAME_LENGTH)';
  i = find( dummy == 0 ,1);
  dummy = char( dummy(1:i-1) );
  varname{v} = dummy;
  clear i dummy
  precisionStr = 'float';
  if precision==8, precisionStr = 'double'; end

  %  Read the variable field depending on what type it is
  if vartype(v) == CV_SCALAR,
        for b=1:nbrBlocks,
          nn=ncvs(b,:);
          tmp = fread(fid, prod(nn), precisionStr);
          tmp = permute( reshape(tmp, nn(3),nn(2),nn(1)) ,[3 2 1]);
          data{v}{b} = tmp;
        end
  elseif vartype(v) == CV_VECTOR,
    for d=1:3,
        for b=1:nbrBlocks,
          nn=ncvs(b,:);
          tmp = fread(fid, prod(nn), precisionStr);
          tmp = permute( reshape(tmp, nn(3),nn(2),nn(1)) ,[3 2 1]);
          data{v}{b}(1:nn(1),1:nn(2),1:nn(3),d) = tmp;
        end
    end
  elseif vartype(v) == CV_TENSOR,
    for d1=1:3,
    for d2=1:3,
        for b=1:nbrBlocks,
          nn=ncvs(b,:);
          tmp = fread(fid, prod(nn), precisionStr);
          tmp = permute( reshape(tmp, nn(3),nn(2),nn(1)) ,[3 2 1]);
          data{v}{b}(1:nn(1),1:nn(2),1:nn(3),d1,d2) = tmp;
        end
    end
    end
  elseif vartype(v) == FA_SCALAR,
      for faceDir=1:3,
        for b=1:nbrBlocks,
          nn=ncvs(b,:); nn(faceDir)=nn(faceDir)+1;
          tmp = fread(fid, prod(nn), precisionStr);
          tmp = permute( reshape(tmp, nn(3),nn(2),nn(1)) ,[3 2 1]);
          if faceDir==1,      data{v}{b}.i = tmp;
          elseif faceDir==2,  data{v}{b}.j = tmp;
          else,               data{v}{b}.k = tmp;
          end
        end
      end
  elseif vartype(v) == FA_VECTOR,
    for d=1:3,
      for faceDir=1:3,
        for b=1:nbrBlocks,
          nn=ncvs(b,:); nn(faceDir)=nn(faceDir)+1;
          tmp = fread(fid, prod(nn), precisionStr);
          tmp = permute( reshape(tmp, nn(3),nn(2),nn(1)) ,[3 2 1]);
          if faceDir==1,      data{v}{b}.i(1:nn(1),1:nn(2),1:nn(3),d) = tmp;
          elseif faceDir==2,  data{v}{b}.j(1:nn(1),1:nn(2),1:nn(3),d) = tmp;
          else,               data{v}{b}.k(1:nn(1),1:nn(2),1:nn(3),d) = tmp;
          end
        end
      end
    end
  elseif vartype(v) == FA_TENSOR,
    for d1=1:3,
    for d2=1:3,
      for faceDir=1:3,
        for b=1:nbrBlocks,
          nn=ncvs(b,:); nn(faceDir)=nn(faceDir)+1;
          tmp = fread(fid, prod(nn), precisionStr);
          tmp = permute( reshape(tmp, nn(3),nn(2),nn(1)) ,[3 2 1]);
          if faceDir==1,      data{v}{b}.i(1:nn(1),1:nn(2),1:nn(3),d1,d2) = tmp;
          elseif faceDir==2,  data{v}{b}.j(1:nn(1),1:nn(2),1:nn(3),d1,d2) = tmp;
          else,               data{v}{b}.k(1:nn(1),1:nn(2),1:nn(3),d1,d2) = tmp;
          end
        end
      end
    end
    end
  elseif vartype(v) == NO_VECTOR,
    for d=1:3,
        for b=1:nbrBlocks,
          nn=ncvs(b,:)+1;
          tmp = fread(fid, prod(nn), precisionStr);
          tmp = permute( reshape(tmp, nn(3),nn(2),nn(1)) ,[3 2 1]);
          data{v}{b}(1:nn(1),1:nn(2),1:nn(3),d) = tmp;
        end
    end
  else,
    disp(['Have not yet implemented reading of vartype == ',int2str(vartype(v)),' ...'])
    return;
  end  

  %  Try reading the next vartype
  tmp = fread(fid, 1, 'integer*4');
end
clear tmp v precision precisionStr nn b d faceDir

fclose(fid);

