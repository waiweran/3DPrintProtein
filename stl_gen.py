import numpy as np
import struct
import datetime
from stl import mesh
import math


class Model3D:
    """Represents a 3D model to make an STL file."""

    def __init__(self):
        self.spheres = list()
        self.cylinders = list()
        self.data = None

    def add_sphere(self, r, pos, granularity, tags):
        """Adds a sphere to the model with the given radius and position (x,y,z) coordinates.
           Granularity determines the number of faces on each ring around the sphere.
           Tags are used for filtering items in model creation."""
        self.spheres.append(Sphere(r, pos, granularity, tags))

    def add_cylinder(self, r, start, end, granularity, tags):
        """Adds a cylinder to the model with the given radius and start, end (x,y,z) coordinates.
           Granularity determines the number of faces on the cylinder.
           Tags are used for filtering items in model creation."""
        self.cylinders.append(Cylinder(r, start, end, granularity, tags))

    def center(self):
        """Centers the model around (0, 0, 0)"""
        pos_sum = np.zeros(3)
        count = 0
        for sphere in self.spheres:
            pos_sum += np.asarray(sphere.pos)
            count += 1
        pos_avg = pos_sum / count
        for sphere in self.spheres:
            sphere.pos = tuple(np.asarray(sphere.pos) - pos_avg)
        for cylinder in self.cylinders:
            cylinder.start = tuple(np.asarray(cylinder.start) - pos_avg)
            cylinder.end = tuple(np.asarray(cylinder.end) - pos_avg)

    def to_stl(self, exclude_tags=None):
        """Generates spheres and cylinders for the STL based on the shapes added to the model.
           Includes only items with all include_tags, and excludes anything with any exclude_tags.
           Empty lists are ignored for both."""
        if exclude_tags is None:
            exclude_tags = []
        data = list()
        for sphere in self.spheres:
            use = True
            for tag in exclude_tags:
                if tag in sphere.tags:
                    use = False
                    break
            if use:
                data.append(_gen_sphere(sphere))
        for cylinder in self.cylinders:
            use = True
            for tag in exclude_tags:
                if tag in cylinder.tags:
                    use = False
                    break
            if use:
                data.append(_gen_cylinder(cylinder))
        if data:
            self.data = _merge_shapes(data)

    def to_binary(self, name):
        """Converts the given STL data to a binary STL string with the given name"""
        output = mesh.Mesh(self.data)
        data = output.data

        # Create the header
        header = '%s (%s) %s %s' % (
            'Online STL Writer',
            '0.1',
            datetime.datetime.now(),
            name,
        )

        # Make it exactly 80 characters
        header = header[:80].ljust(80, ' ')
        header = bytes(header, 'ascii', 'replace')

        packed = struct.pack('<i', data.size)
        content = data.tobytes()

        return header + packed + content

    def save(self, filename):
        """Writes the given STL data to a file with the given name"""
        output = mesh.Mesh(self.data)
        output.save(filename)


def _merge_shapes(data):
    """Takes a list of STL objects and combines them into one STL object"""
    return np.concatenate(data)


def _gen_cylinder(cylinder):
    """Generates a cylinder with the given radius and start, end (x,y,z) coordinates.
       Granularity determines the number of faces on the cylinder."""
    vector = np.array(cylinder.end) - np.array(cylinder.start)
    axis = vector / np.linalg.norm(vector)
    rad_unsc = np.cross(vector, np.array([1, 1, 0]))
    if rad_unsc[0] == 0 and rad_unsc[1] == 0 and rad_unsc[2] == 0:
        rad_unsc = np.cross(vector, np.array([1, 0, 0]))
    radius = rad_unsc/np.linalg.norm(rad_unsc) * cylinder.radius

    # Points on Circle
    count = 0
    circle_pts = np.zeros((cylinder.granularity, 3))
    for theta in np.linspace(0, 2*math.pi*(1 - 1/cylinder.granularity), cylinder.granularity):
        circle_pts[count, :] = radius*np.cos(theta) + np.cross(axis, radius)*np.sin(theta) + \
                              axis*np.dot(axis, radius)*(1 - np.cos(theta))
        count += 1

    # Cylinder Sides
    count = 0
    coords = np.zeros((cylinder.granularity*2, 3, 3))
    for i in range(0, cylinder.granularity - 1):
        coords[count, :] = np.array([circle_pts[i, :] + cylinder.end,
                                     circle_pts[i, :] + cylinder.start,
                                     circle_pts[i+1, :] + cylinder.start])
        coords[count+1, :] = np.array([circle_pts[i, :] + cylinder.end,
                                      circle_pts[i+1, :] + cylinder.start,
                                      circle_pts[i+1, :] + cylinder.end])
        count += 2
    coords[count, :] = np.array([circle_pts[cylinder.granularity-1, :] + cylinder.end,
                                 circle_pts[cylinder.granularity-1, :] + cylinder.start,
                                 circle_pts[0, :] + cylinder.start])
    coords[count+1, :] = np.array([circle_pts[cylinder.granularity-1, :] + cylinder.end,
                                  circle_pts[0, :] + cylinder.start,
                                  circle_pts[0, :] + cylinder.end])
    count += 2

    # Convert to STL Generator
    output = np.zeros(count, dtype=mesh.Mesh.dtype)
    output['vectors'][:] = coords[:]

    return output


def _gen_sphere(sphere):
    """Generates a sphere with the given radius and position (x,y,z) coordinates.
       Granularity determines the number of faces on each ring around the sphere."""
    coords = np.zeros((sphere.granularity*2 + 2*sphere.granularity*(sphere.granularity-3), 3, 3))
    thetas = np.linspace(0, math.pi, sphere.granularity)
    phis = np.linspace(0, 2*math.pi, sphere.granularity+1)
    count = 0

    # Bottom Triangles
    for j in range(0, len(phis) - 1):
        coords[count, :] = np.array([[sphere.radius, thetas[0], 0],
                                     [sphere.radius, thetas[1], phis[j]],
                                     [sphere.radius, thetas[1], phis[j+1]]])
        count += 1

    # Central Tiles
    for i in range(1, len(thetas) - 2):
        for j in range(0, len(phis) - 1):
            coords[count, :] = np.array([[sphere.radius, thetas[i], phis[j]],
                                         [sphere.radius, thetas[i+1], phis[j]],
                                         [sphere.radius, thetas[i], phis[j+1]]])
            coords[count+1, :] = np.array([[sphere.radius, thetas[i+1], phis[j+1]],
                                          [sphere.radius, thetas[i], phis[j+1]],
                                          [sphere.radius, thetas[i+1], phis[j]]])
            count += 2

    # Top Triangles
    for j in range(0, len(phis) - 1):
        coords[count, :] = np.array([[sphere.radius, thetas[-2], phis[j+1]],
                                     [sphere.radius, thetas[-2], phis[j]],
                                     [sphere.radius, thetas[-1], 0]])
        count += 1

    # Map to Cartesian
    cartesian = np.zeros((count, 3, 3))
    cartesian[:, :, 0] = sphere.radius * np.sin(coords[:, :, 1]) * np.cos(coords[:, :, 2])
    cartesian[:, :, 1] = sphere.radius * np.sin(coords[:, :, 1]) * np.sin(coords[:, :, 2])
    cartesian[:, :, 2] = sphere.radius * np.cos(coords[:, :, 1])

    # Convert to STL Generator
    output = np.zeros(count, dtype=mesh.Mesh.dtype)
    output['vectors'][:] = cartesian[:]

    output['vectors'][:, :, 0] += sphere.pos[0]
    output['vectors'][:, :, 1] += sphere.pos[1]
    output['vectors'][:, :, 2] += sphere.pos[2]

    return output


class Sphere:
    """Represents a sphere"""

    def __init__(self, radius, pos, granularity, tags):
        self.radius = radius
        self.pos = pos
        self.granularity = granularity
        self.tags = tags


class Cylinder:
    """Represents a cylinder"""

    def __init__(self, radius, start, end, granularity, tags):
        self.radius = radius
        self.start = start
        self.end = end
        self.granularity = granularity
        self.tags = tags
