save_path = "/media/sf_Virtual_Drive/Winter 2021/Wake Visualization/";
vtk_save_name="test_me";

# Define the center of the grid.
center = [4.0,0.0,0.0];
x_length = 2.0;
y_length = 4.0;
z_length = 2.0;

r = center[1];
v = 10;

t_total = 1.0;
n_steps = 200;



for i in 1:n_steps;

    t = i*(t_total/n_steps) # This is the real time. 

    x1=r*cos(t*v/r)-x_length/2;
    x2=r*cos(t*v/r)+x_length/2;
    y1=r*sin(t*v/r)+y_length/2;
    y2=r*sin(t*v/r)-y_length/2;

    circle = [[min(x1,x2),min(y1,y2),-1],[max(x1,x2),max(y1,y2),1]]

    fdom = gt.Grid(circle[1],circle[2],convert(Array{Int64,1},[1,1,1]))

    gt.save(fdom,"$save_path$vtk_save_name";num=i)

end