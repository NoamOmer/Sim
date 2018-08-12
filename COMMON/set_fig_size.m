function set_fig_size(fig,paper_size)
if nargin==1
    paper_size=fig;
    fig=gcf;
end
root_units=get(0,'units');
set(0,'units','inch')
screen_size=get(0,'screensize');
set(0,'units',root_units)
figure_pos=paper_size+...
    [(screen_size(3)-paper_size(3))/2 (screen_size(4)-paper_size(4))/2 0 0]
set(fig,'units','inch','pos',figure_pos,...
    'paperunits','inch','paperpos',paper_size,'color','w')