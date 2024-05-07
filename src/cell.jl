module cellModule

export Cell

struct Cell
    """
    the base type for other cells: vertices, edges, faces, ...
    """
    dimension::Int

    function Cell(dimension)
        new(dimension)
    end
end
end
