type Triangle
  p1::Vertex
  p2::Vertex
  p3::Vertex
  #TODO order the points counter clockwise. Maybe have a sort function for Vertex
  Triangle(p1,p2,p3) = new(p1,p2,p3)
end
