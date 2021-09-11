extern crate cairo;
extern crate rand;

fn draw_bounds(context: & cairo::Context, points: &[(f64, f64)], r:f64, g:f64, b:f64)
{
    let mut min_x = points[0].0;
    let mut min_y = points[0].1;
    let mut max_x = points[0].0;
    let mut max_y = points[0].1;
    
    for (x, y) in points
    {
        if *x < min_x {min_x = *x};
        if *y < min_y {min_y = *y};
        if *x > max_x {max_x = *x};
        if *y > max_y {max_y = *y};

    }
    
    context.rectangle(min_x as f64, min_y as f64, (max_x - min_x) as f64, (max_y - min_y) as f64);
    context.set_source_rgb(r, g, b);
    context.set_line_width(0.01);
    context.stroke().unwrap();
}

fn main() {
    let surface = cairo::PdfSurface::new(100.0, 100.0, "out.pdf").unwrap();
    let context = cairo::Context::new(&surface).unwrap();
    
    use rand::Rng;
    let mut rng = rand::XorShiftRng::new_unseeded();
    
    let mut points = Vec::new();
    const NUM_POINTS : u32 = 10000;
    const NUM_CHILD_NODES : u32 = 8;
    
    let tree_height = (NUM_POINTS as f64).log(NUM_CHILD_NODES as f64).ceil() as i32;
    let nodes_on_last_interior_level = (NUM_CHILD_NODES).checked_pow((tree_height - 1) as u32).unwrap();
    
    println!("target tree height: {}", tree_height);
    println!("num nodes on last interior level: {}", nodes_on_last_interior_level);
    let max_num_children_per_node = (NUM_POINTS as f64 / nodes_on_last_interior_level as f64).ceil() as u32;
    println!("num children per node: {}", max_num_children_per_node);
    // NUM_POINTS = nodes_will_full_child_count * max_num_children_per_node + (nodes_on_last_interior_level - nodes_will_full_child_count) * (max_num_children_per_node-1)
    let nodes_with_full_child_count = NUM_POINTS - nodes_on_last_interior_level * (max_num_children_per_node-1);
    println!("{} = {} * {} + {} * {}", NUM_POINTS, nodes_with_full_child_count, max_num_children_per_node, nodes_on_last_interior_level - nodes_with_full_child_count, max_num_children_per_node-1);
    
    
    for _i in 0..NUM_POINTS
    {
        points.push((rng.next_f64()*100.0, rng.next_f64()*100.0));
    }
    
    points.sort_by(|lhs, rhs|{ lhs.0.partial_cmp(&rhs.0).unwrap().then_with(|| lhs.1.partial_cmp(&rhs.1).unwrap())});
 
 
    for i in 0..10
    {
        let mut slice = points[i*100..std::cmp::min((i+1)*100, 1000)].to_vec();
        draw_bounds(& context, &slice, 0.1, 0.1, (i as f64) *0.1);
        
        slice.sort_by(|lhs, rhs|{ lhs.1.partial_cmp(&rhs.1).unwrap()});
        for i in 0..10
        {
            let the_box = &slice[i*10..(i+1)*10];
            draw_bounds( &context, the_box, 0.1, (i as f64) *0.1, 0.1);
        }
    }
//    let mut slices = Vec::new();
    let mut r:f64 = 0.0;
    for (x,y) in points
    {
        context.rectangle(x-0.05, y-0.05, 0.1, 0.1);
        context.set_line_width(0.1);
        context.set_source_rgb(r, 0.0, 0.0);
        r+= 0.001;
        context.fill().unwrap();
    }
//    context.move_to(0.0, 0.0);
//    context.line_to(512.0, 512.0);
    println!("Hello, world!");
}
