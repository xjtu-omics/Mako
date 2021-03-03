
class Segment():
    # x indicates read, while y ref
    def __init__(self, x_start, y_start, length, forward, seg_id):
        self.x_start = x_start
        self.y_start = y_start
        self.seg_length = length
        self.seg_forward = forward
        self.seg_id = seg_id

        self.seg_flag = ""
        if forward:
            self.x_end = self.x_start + (self.seg_length - 1)
        else:
            self.x_end = self.x_start - (self.seg_length - 1)
        self.y_end = self.y_start + (length - 1)

    def set_flag(self, flag):
        self.seg_flag = flag

    def set_x_end(self, xEnd):
        self.x_end = xEnd
    
    def set_y_end(self, yEnd):
        self.y_end = yEnd
    
    def set_length(self, length):
        self.seg_length = length
    
    def set_x_start(self, xStart):
        self.x_start = xStart
    
    def set_y_start(self, yStart):
        self.y_start = yStart
    
    def set_forward(self, forward):
        self.seg_forward = forward
    
    def set_seg_id(self, segId):
        self.seg_id = segId

    def to_string(self):
        coord_x = "[" + str(self.x_start) + "," + str(self.x_end) + "]"
        coord_y = "[" + str(self.y_start) + "," + str(self.y_end) + "]"
        return str(self.seg_id) + "-" + coord_x + "-" + coord_y + "-" + str(self.seg_forward) + "-" + str(self.seg_length) + "-" + str(self.seg_flag)

