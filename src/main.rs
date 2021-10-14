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

fn sort_key( value: &(f64,f64), dimension: usize) -> f64 {
    return if dimension % 2 == 0
        { value.0 }
    else
        { value.1 }
}

fn main() {
    let surface = cairo::PdfSurface::new(100.0, 100.0, "out.pdf").unwrap();
    let context = cairo::Context::new(&surface).unwrap();
    
    use rand::Rng;
    let mut rng = rand::XorShiftRng::new_unseeded();
    
    let mut points = Vec::new();
    const NUM_POINTS : usize = 10000;
    const NUM_CHILD_NODES : usize = 18;
    
    let tree_height = (NUM_POINTS as f64).log(NUM_CHILD_NODES as f64).ceil() as i32;
    println!("target tree height: {}", tree_height);

    let max_num_children_per_node = (NUM_POINTS as f64).powf(1.0/tree_height as f64).ceil();
    println!("num children per node: {}", max_num_children_per_node);
    
    for _i in 0..NUM_POINTS
    {
        points.push((rng.next_f64()*100.0, rng.next_f64()*100.0));
    }
    
    let num_blocks = NUM_POINTS as f64 / max_num_children_per_node as f64;
    println!("Total num blocks: {}", num_blocks);
    let num_nodes_per_slice = NUM_POINTS as f64 / (num_blocks.sqrt());
    let num_slices = (NUM_POINTS as f64 / num_nodes_per_slice).ceil() as usize ;
    println!("#slices: {}", num_slices);
 
    let num_nodes_per_slice = (NUM_POINTS as f64 / num_slices as f64).ceil() as usize;
    println!("Nodes per slice: {}", num_nodes_per_slice);

    let mut start_index : usize = 0;
    let mut num_blocks_created = 0;

    points.sort_by(|lhs, rhs|{ sort_key(lhs, 0).partial_cmp(&sort_key(rhs,0)).unwrap()});

    for i in 0..num_slices
    {
       
        let nodes_in_this_slice = if (num_slices - i) * (num_nodes_per_slice-1) < (NUM_POINTS - start_index) { num_nodes_per_slice } else {num_nodes_per_slice - 1};
        
        let mut slice = points[start_index..start_index + nodes_in_this_slice].to_vec();
        draw_bounds(& context, &slice, 0.1, 0.1, (i as f64) *0.1);
        
        slice.sort_by(|lhs, rhs|{ sort_key(lhs, 1).partial_cmp(&sort_key(rhs, 1)).unwrap()});


        let num_blocks = ((slice.len() as f64) / (max_num_children_per_node as f64)).ceil() as usize;
        let num_nodes_per_block = ((slice.len() as f64) / num_blocks as f64).ceil() as usize;
        println!("slice {}, {} nodes, {} blocks, {} nodes per block", i, slice.len(), num_blocks, num_nodes_per_block);
        
        let mut block_start_index : usize = 0;

        for i in 0..num_blocks
        {
            let nodes_in_this_block = if (num_blocks - i) * (num_nodes_per_block-1) < (slice.len() - block_start_index) { num_nodes_per_block } else {num_nodes_per_block - 1};

            let the_box = &slice[block_start_index..block_start_index + nodes_in_this_block];
            draw_bounds( &context, the_box, 0.1, (i as f64) *0.1, 0.1);
            block_start_index += nodes_in_this_block;
            num_blocks_created+=1;
        }
        start_index += nodes_in_this_slice;
    }

    let mut r:f64 = 0.0;
    for (x,y) in points
    {
        context.rectangle(x-0.05, y-0.05, 0.1, 0.1);
        context.set_line_width(0.1);
        context.set_source_rgb(r, 0.0, 0.0);
        r+= 0.001;
        context.fill().unwrap();
    }
    println!("Created a total of {} blocks", num_blocks_created);
}
