function [x_cv, x_no, x_fa, ncvs, ...
          bndNames, bnd_block_begin_end, ...
          int_block_begin_end_A, int_block_begin_end_B, int_offset, int_angle] = readTurtleGrid(filename, skipCoords)
%
%  Read grid in Turtle format.
%  Store coordinates as (i,j,k,dir).
%  For faces, store as x_fa.i etc.
%
%  If called with skipCoords, then those fields will be empty on return (but the function is then much faster)
%

%  Parameters from Turtle_parameters.h
MAGICNUMBER    = 751123;
IO_NAME_LENGTH =    200;


if exist('skipCoords', 'var')==0,
  skipCoords = 0;
end


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


%  Read the header
version = fread(fid, 1, 'integer*4');
if version < 10000*(3) + 100*(1) + (1),
  fclose(fid);
  %disp('Error: this version only reads Tortuga files from version 3.1 and on')
  disp('Re-directing to the old version of readTurtleGrid.m')
  [x_cv, x_no, x_fa, ncvs, ...
    bndNames, bnd_block_IJK, ...
    int_block_IJK_A, int_block_IJK_B, int_offset] = readTurtleGrid_old(filename);
    bnd_block_begin_end = [];
    int_block_begin_end_A = [];
    int_block_begin_end_B = [];
    clear bnd_block_IJK int_block_IJK_A int_block_IJK_B
  return
end
nbrBlocks = fread(fid, 1, 'integer*4');
nbrInterfaces = fread(fid, 1, 'integer*4');
nbrBoundaries = fread(fid, 1, 'integer*4');

%  Read the size of each block
for b=1:nbrBlocks,
  ncvs(b,:) = fread(fid, 3, 'integer*4');
end

%  Read interface boundaries
int_block_begin_end_A=[]; int_block_begin_end_B=[]; int_offset=[]; int_angle=[];
for bnd=1:nbrInterfaces,
  int_block_begin_end_A{bnd} = fread(fid, 7, 'integer*4')';
  int_block_begin_end_B{bnd} = fread(fid, 7, 'integer*4')';
  int_offset{bnd} = fread(fid, 3, 'double')';
  if version < 10000*(3) + 100*(1) + (47),
    int_angle{bnd} = 0;
    if bnd==1,
      disp('Setting interface angle to zero since reading old format...')
    end
  else,
    int_angle{bnd} = fread(fid, 1, 'double');
  end
end

%  Regular boundaries
bndNames=[]; bnd_block_begin_end=[];
for bnd=1:nbrBoundaries,
  tmp = fread(fid, IO_NAME_LENGTH)';
  i = find( tmp <= 32 ,1);
  bndNames{bnd} = char( tmp(1:i-1) );
  nbrBlockFaces = fread(fid, 1, 'integer*4');
  for i=1:nbrBlockFaces,
    bnd_block_begin_end{bnd}(i,:) = fread(fid, 7, 'integer*4');
  end
end

%  Possibly end here
if skipCoords > 0,
  disp('  skipping the coordinates...')
  x_cv = [];
  x_fa = [];
  x_no = [];
  fclose(fid);
  return;
end


%  Node coordinates
%  Stored in kji-ordering, block by block
for d=1:3,
  for b=1:nbrBlocks, for i=1:ncvs(b,1)+1, for j=1:ncvs(b,2)+1,
    tmp = fread(fid, ncvs(b,3)+1, 'double');
    x_no{b}(i,j,:,d) = permute(tmp, [3 2 1]);
  end;end;end
end

fclose(fid);



%  Compute cv centroids
for b=1:nbrBlocks,
  i=1:ncvs(b,1); j=1:ncvs(b,2); k=1:ncvs(b,3);
  x_cv{b} = 0.125*( x_no{b}(i  ,j  ,k  ,:) + x_no{b}(i+1,j  ,k  ,:) + x_no{b}(i  ,j+1,k  ,:) + x_no{b}(i+1,j+1,k  ,:) ...
                  + x_no{b}(i  ,j  ,k+1,:) + x_no{b}(i+1,j  ,k+1,:) + x_no{b}(i  ,j+1,k+1,:) + x_no{b}(i+1,j+1,k+1,:) );
end

%  Compute face centroids
for b=1:nbrBlocks,
  i=1:ncvs(b,1); j=1:ncvs(b,2); k=1:ncvs(b,3);
  I=1:ncvs(b,1)+1; J=1:ncvs(b,2)+1; K=1:ncvs(b,3)+1;
  x_fa{b}.i = 0.25*( x_no{b}(I  ,j  ,k  ,:) + x_no{b}(I  ,j+1,k  ,:) + x_no{b}(I  ,j  ,k+1,:) + x_no{b}(I  ,j+1,k+1,:) );
  x_fa{b}.j = 0.25*( x_no{b}(i  ,J  ,k  ,:) + x_no{b}(i+1,J  ,k  ,:) + x_no{b}(i  ,J  ,k+1,:) + x_no{b}(i+1,J  ,k+1,:) );
  x_fa{b}.k = 0.25*( x_no{b}(i  ,j  ,K  ,:) + x_no{b}(i+1,j  ,K  ,:) + x_no{b}(i  ,j+1,K  ,:) + x_no{b}(i+1,j+1,K  ,:) );
end

