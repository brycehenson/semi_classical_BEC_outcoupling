function hist_2d_data(data_in,opts)
% plot a histogram 'slice' or integeration of xyz data
% opts.cords_to_plot
% opts.cords_labels
% opts.blur_size
% opts.hist_bounds

% optional
% opts.filter_in_int_axis
% opts.bin_factor;
% opts.std_norm_bounds
    % opts.mean_cen;
%opts.cords_units={'m','m','m'};
%opts.scale_fac=[1,1,1]*1e-6;
%opts.mean_marker
% opts.cmap
% opts.norm_den_to_num_in
% opts.scale_den_disp_fac
% opts.save_hist_as_image
% opts.save_hist_name
% opts.norm_den_unity



cords_to_plot=opts.cords_to_plot;
cords_labels=opts.cords_labels;

if ~isfield(opts,'filter_in_other_axis')
    opts.filter_in_other_axis=false;
end
if ~isfield(opts,'std_norm_bounds')
    opts.std_norm_bounds=false;
end
if ~isfield(opts,'opts.bin_factor')
    opts.bin_factor=5;
end
if ~isfield(opts,'disp_scale')
    opts.disp_scale=[1,1,1]; 
end
if ~isfield(opts,'mean_marker')
    opts.mean_marker=false;
end
if ~isfield(opts,'norm_den_to_num_in')
    opts.norm_den_to_num_in=false;
end
if ~isfield(opts,'scale_den_disp_fac')
    opts.scale_den_disp_fac=false;
end
if ~isfield(opts,'norm_den_unity')
    opts.norm_den_unity=false;
end
if ~isfield(opts,'save_hist_as_image')
    opts.save_hist_as_image=false;
end

% check that size of data is right
if max(cords_to_plot)>size(data_in,1)
    error('max cord to plot exceeds data size')
end
% remove nan from data
data_in=data_in(~any(isnan(data_in),2),:);
data_mean=mean(data_in,1)';
num_counts=size(data_in,1);

if opts.std_norm_bounds
    if ~isfield(opts,'mean_cen')
        opts.mean_cen=true;
    end
    % define the bounds in terms of the standard deviation of the data in each axis
    data_std=std(data_in,1)';
    
    bound_size_in_sd=col_vec(opts.hist_bounds);
    hist_bounds=[[-1,1];[-1,1];[-1,1]].*repmat(data_std.*(bound_size_in_sd),[1,2]);
    if  opts.mean_cen
        hist_bounds=hist_bounds+repmat(data_mean,[1,2]);
    end
    blur_size=data_std(cords_to_plot)'.*opts.blur_size(cords_to_plot);
else 
    hist_bounds=opts.hist_bounds;
    blur_size=opts.blur_size(cords_to_plot);
end

hist_bins=ones(3,1);
bounds_range_tmp=diff(hist_bounds,1,2);
bounds_range_tmp=bounds_range_tmp(cords_to_plot);
hist_bins(cords_to_plot)=ceil((bounds_range_tmp./blur_size')*opts.bin_factor);


stfig('2d hist'); %'add_stack',1
%[[xmin,xmax];[ymin,ymax];[zmin,zmax]]

hist_edges=cell(1,3);
for ii=1:3
    hist_edges{ii}=linspace(hist_bounds(ii,1),hist_bounds(ii,2),hist_bins(ii));
end
lin_bin_size=range(hist_bounds,2)./hist_bins;

if opts.filter_in_other_axis
    all_cord=1:3;
    int_cord=all_cord;
    int_cord(cords_to_plot)=[];
    mask_int_axis=data_in(:,int_cord)>hist_bounds(int_cord,1)  &  data_in(:,int_cord)<hist_bounds(int_cord,2);
    data_in=data_in(mask_int_axis,:);
end

bin_area=prod(lin_bin_size(cords_to_plot));
[counts,centers]=hist3(data_in(:,cords_to_plot),'edges',{hist_edges{cords_to_plot(1)},hist_edges{cords_to_plot(2)}});
% sum(sum(counts==0))
% max(max(counts))
counts=counts/bin_area; 

if opts.norm_den_to_num_in
    counts=counts/num_counts;
end

if opts.scale_den_disp_fac
    % normalize so the density is in the same units
    den_disp_scale=prod(opts.disp_scale(cords_to_plot));
    counts=counts*den_disp_scale;
end

if  any(blur_size~=0)
    kernel_sigma_in_bins=blur_size./lin_bin_size(cords_to_plot);
    kernel_size=ceil(kernel_sigma_in_bins*5); %size of kernel
    blur_kernel=custom_gaussian_function_2d(kernel_size, kernel_sigma_in_bins, [], [], 0, 1);
    counts=imfilter(counts, blur_kernel, 'same');
end

if opts.norm_den_unity
    counts=counts/max(counts(:));
end

%imagesc seems to plot the wrong way round so we transpose here
imagesc(centers{1}/opts.disp_scale(cords_to_plot(1)),...
    centers{2}/opts.disp_scale(cords_to_plot(2)),counts')

if ~isfield(opts,'cmap')
    opts.cmap=viridis;
end
colormap(opts.cmap)
set(gca,'Ydir','normal')
set(gcf,'Color',[1 1 1]);

if isfield(opts,'cbar_lims')
    caxis(opts.cbar_lims);
end

if opts.mean_marker
hold on
plot(data_mean(cords_to_plot(1))/opts.disp_scale(cords_to_plot(1)),...
    data_mean(cords_to_plot(2))/opts.disp_scale(cords_to_plot(2)),...
    'rx','MarkerSize',12,'LineWidth',2)
hold off
end


prefixes={'','',''};
% conver the opts.disp_scale to a unit prefix
prefix_table={1e15,'P';1e12,'T';1e9,'G';1e6,'M';1e3,'k';1e-3,'m';1e-6,'$\mu$';1e-9,'n';1e-12,'p';1e-15,'f'};
factors=[prefix_table{:,1}];
factors_eps=eps(factors)*10;
for ii=1:3
    fac_diff=abs(opts.disp_scale(ii)-factors);
    fac_diff=fac_diff<factors_eps;
    idx=find(fac_diff);
    if numel(idx)>1
        error('should not get multiple index')
    end
    if numel(idx)==1
        prefixes{ii}= prefix_table{idx,2};
    end
end



xlabel(sprintf('%s (%s%s)',...
    opts.cords_labels{cords_to_plot(1)},...
    prefixes{cords_to_plot(1)},...
    opts.cords_units{cords_to_plot(1)}))
ylabel(sprintf('%s (%s%s)',...
    opts.cords_labels{cords_to_plot(2)},...
    prefixes{cords_to_plot(2)},...
    opts.cords_units{cords_to_plot(2)}))

c=colorbar;
set(c,'LineWidth',1.1,'TickLength',[0.025]);
if opts.norm_den_unity
    den_units='arb. u.';
elseif opts.scale_den_disp_fac
    if numel(unique(opts.disp_scale(cords_to_plot)))==1
        den_units=sprintf('(%s%s)$^{-2}$',...
                            prefixes{cords_to_plot(1)},...
                            opts.cords_units{cords_to_plot(1)});
    else
        den_units=sprintf('(%s%s%s%s)$^{-1}$',...
                            prefixes{cords_to_plot(1)},...
                            opts.cords_units{cords_to_plot(1)},...
                            prefixes{cords_to_plot(2)},...
                            opts.cords_units{cords_to_plot(2)});
    end
else
    if strcmp(opts.cords_units{cords_to_plot(1)},opts.cords_units{cords_to_plot(2)})
        den_units=sprintf('(%s)$^{-2}$',...
                            opts.cords_units{cords_to_plot(1)});
    else
        den_units=sprintf('(%s%s)$^{-1}$',...
                            opts.cords_units{cords_to_plot(1)},...
                            opts.cords_units{cords_to_plot(2)});
    end
end
%den_units=strrep(den_units,'$$','');
c.Label.Interpreter = 'latex';
c.Label.String=sprintf('Density (%s)',den_units);
daspect([1,1,1])

if opts.save_hist_as_image
    if ~isfield(opts,'save_hist_name')
        opts.save_hist_name='out_hist.tiff';
    end
    unit16max=2^16-1;
    norm_counts=uint16(counts/max(max(counts))*unit16max);
    norm_counts=transpose(norm_counts);
    norm_counts=flipud(norm_counts);
    cmap_unit16=interp_colormap(opts.cmap,linspace(0,1,unit16max));
    rgbcounts= ind2rgb(norm_counts, cmap_unit16);
    imwrite(rgbcounts,opts.save_hist_name)
end

end