"""Module with objects that could be read or written to netCDF"""
import copy
import os
from dataclasses import dataclass
from typing import ClassVar, List, Optional, Tuple

from netCDF4 import Dataset, Variable  # type: ignore
from numpy import ma
import numpy.typing as npt

from palm_csd.csd_config import CSDConfig, CSDConfigAttributes, CSDConfigDomain, CSDConfigInput


@dataclass
class NCDFDimension:
    """A dimension that could written to netCDF. A corresponding dimension variable is included.
    Its values can be stored in the attribute values or supplied when written.
    """

    name: str
    datatype: str

    long_name: Optional[str] = None
    units: Optional[str] = None
    standard_name: Optional[str] = None
    values: Optional[npt.NDArray] = None

    _metadata: ClassVar[List[str]] = ["long_name", "units", "standard_name"]

    @property
    def size(self) -> int:
        """Get the size of the dimension. It is derived from the values attribute."""
        if self.values is None:
            raise ValueError("Values of dimension " + self.name + " not defined.")
        return self.values.size

    def __len__(self) -> int:
        return self.size

    def to_dataset(self, nc_data: Dataset, values: Optional[ma.MaskedArray] = None) -> None:
        """Add dimension and its dimension variable to a netCDF Dataset if it does not already
        include it. The values are either supplied in the function call or taken from the values
        attribute.
        """
        if self.name not in nc_data.dimensions:
            if values is not None:
                to_write = values
            elif self.values is not None:
                to_write = self.values
            else:
                raise ValueError("Values of dimension " + self.name + " not defined.")

            print("Writing dimension " + self.name + " to file...")

            nc_data.createDimension(self.name, len(to_write))

            nc_var = nc_data.createVariable(self.name, self.datatype, self.name)
            nc_var[:] = to_write

            for attr in self._metadata:
                attr_value = getattr(self, attr)
                if attr_value is not None:
                    nc_var.setncattr(attr, attr_value)


@dataclass
class NCDFVariable:
    """A variable that could written to netCDF. It includes its metadata and dimensions.
    Its values can be stored in the attribute values or supplied when written. A default filename
    can be also stored.
    """

    name: str
    dimensions: Tuple[NCDFDimension, ...]
    datatype: str
    fillvalue: float
    long_name: str
    units: str

    values: Optional[npt.NDArray] = None

    coordinates: Optional[str] = None
    grid_mapping: Optional[str] = None
    lod: Optional[int] = None
    res_orig: Optional[float] = None
    standard_name: Optional[str] = None

    filename: Optional[str] = None

    _metadata: ClassVar[List[str]] = [
        "long_name",
        "units",
        "standard_name",
        "res_orig",
        "lod",
        "coordinates",
        "grid_mapping",
    ]

    def to_file(
        self, values: Optional[ma.MaskedArray] = None, filename: Optional[str] = None
    ) -> None:
        """Write variable to a netCDF file, which is openend and closed. The file is either
        specified in the function call or taken from the default filename. If the variable was not
        yet added the file, it is added with its dimensions, otherwise its values are overwritten.
        The values are either supplied in the function call or taken from the values attribute.
        """

        if values is not None:
            to_write = values
        elif self.values is not None:
            to_write = self.values
        else:
            raise ValueError("Values of variable " + self.name + " not defined.")

        if filename is not None:
            to_file = filename
        elif self.filename is not None:
            to_file = self.filename
        else:
            raise ValueError("Output filename for variable " + self.name + " not defined.")

        try:
            nc_data = Dataset(to_file, "a", format="NETCDF4")
        except FileNotFoundError:
            print("Error. Could not open file: " + to_file + ". Aborting...")
            raise

        print("Writing array " + self.name + " to file...")
        if self.name not in nc_data.variables:
            for nc_dim in self.dimensions:
                nc_dim.to_dataset(nc_data)

            nc_var = nc_data.createVariable(
                self.name,
                self.datatype,
                (o.name for o in self.dimensions),
                fill_value=self.fillvalue,
            )

            for attr in self._metadata:
                attr_value = getattr(self, attr)
                if attr_value is not None:
                    nc_var.setncattr(attr, attr_value)

        else:
            nc_var = nc_data.variables[self.name]

        if len(self.dimensions) == 1:
            nc_var[:] = to_write
        elif len(self.dimensions) == 2:
            nc_var[:, :] = to_write
        elif len(self.dimensions) == 3:
            nc_var[:, :, :] = to_write
        else:
            raise NotImplementedError

        nc_data.close()

    def from_file(self, filename: Optional[str] = None) -> ma.MaskedArray:
        """Return the variable from a netCDF file. The file is either specified in the
        function call or taken from the default filename.
        """

        if filename is not None:
            from_file = filename
        elif self.filename is not None:
            from_file = self.filename
        else:
            raise ValueError("Input filename for variable " + self.name + " not defined.")

        try:
            nc_data = Dataset(from_file, "r", format="NETCDF4")
        except FileNotFoundError:
            print("Error. Could not open file: " + from_file + ". Aborting...")
            raise

        tmp_array = nc_data.variables[self.name][:]
        nc_data.close()

        return tmp_array


@dataclass
class NCDFCoordinateReferenceSystem:
    """A coordinate reference system that can be written to a netCDF file."""

    long_name: str
    grid_mapping_name: str
    semi_major_axis: float
    inverse_flattening: float
    longitude_of_prime_meridian: float
    longitude_of_central_meridian: float
    scale_factor_at_central_meridian: float
    latitude_of_projection_origin: float
    false_easting: float
    false_northing: float
    spatial_ref: str
    units: str
    epsg_code: str

    filename: Optional[str] = None

    def to_file(self, filename: Optional[str] = None) -> None:
        """Write CRS to a netCDF file, which is openend and closed. The file is either specified
        in the function call or taken from the default filename.
        """

        if filename is not None:
            to_file = filename
        elif self.filename is not None:
            to_file = self.filename
        else:
            raise ValueError("Output filename for CRS not defined.")

        try:
            nc_data = Dataset(to_file, "a", format="NETCDF4")
        except FileNotFoundError:
            print("Error. Could not open file: " + to_file + ". Aborting...")
            raise

        print("Writing crs to file...")

        nc_var = nc_data.createVariable("crs", "i")

        nc_var.long_name = self.long_name
        nc_var.grid_mapping_name = self.grid_mapping_name
        nc_var.semi_major_axis = self.semi_major_axis
        nc_var.inverse_flattening = self.inverse_flattening
        nc_var.longitude_of_prime_meridian = self.longitude_of_prime_meridian
        nc_var.longitude_of_central_meridian = self.longitude_of_central_meridian
        nc_var.scale_factor_at_central_meridian = self.scale_factor_at_central_meridian
        nc_var.latitude_of_projection_origin = self.latitude_of_projection_origin
        nc_var.false_easting = self.false_easting
        nc_var.false_northing = self.false_northing
        nc_var.spatial_ref = self.spatial_ref
        nc_var.units = self.units
        nc_var.epsg_code = self.epsg_code

        nc_data.close()


class NCDFDomain:
    """A domain that stores all its configurations, output dimensions and variables."""

    name: str

    config: CSDConfigDomain
    input_config: CSDConfigInput
    attributes: CSDConfigAttributes

    # TODO: Python 3.11: Use Self to indicate same type as class
    parent: Optional["NCDFDomain"]

    rotation_angle: float

    x: NCDFDimension
    y: NCDFDimension
    z: NCDFDimension
    zlad: NCDFDimension

    nsurface_fraction: NCDFDimension
    nbuilding_pars: NCDFDimension
    nvegetation_pars: NCDFDimension
    nwater_pars: NCDFDimension

    lat: NCDFVariable
    lon: NCDFVariable

    x_UTM: NCDFVariable
    y_UTM: NCDFVariable

    E_UTM: NCDFVariable
    N_UTM: NCDFVariable

    zt: NCDFVariable

    buildings_2d: NCDFVariable
    building_id: NCDFVariable
    building_type: NCDFVariable
    buildings_3d: NCDFVariable

    surface_fraction: NCDFVariable

    vegetation_type: NCDFVariable
    pavement_type: NCDFVariable
    water_type: NCDFVariable
    soil_type: NCDFVariable
    street_type: NCDFVariable
    street_crossing: NCDFVariable

    building_pars: NCDFVariable
    vegetation_pars: NCDFVariable
    water_pars: NCDFVariable

    lad: NCDFVariable
    bad: NCDFVariable
    tree_id: NCDFVariable
    tree_type: NCDFVariable

    crs: NCDFCoordinateReferenceSystem

    def __init__(
        self,
        name: str,
        config: CSDConfig,
        parent: Optional["NCDFDomain"] = None,
    ) -> None:
        self.name = name

        # configurations
        self.config = config.domain_dict[name]
        self.input_config = config.input_of_domain(self.config)

        self.attributes = copy.copy(config.attributes)  # copy for domain-dependent values

        self.parent = parent

        self.rotation_angle = config.settings.rotation_angle

        # dimensions
        self.x = NCDFDimension(
            name="x",
            datatype="f4",
            standard_name="projection_x_coordinate",
            long_name="x",
            units="m",
        )
        self.y = NCDFDimension(
            name="y",
            datatype="f4",
            standard_name="projection_y_coordinate",
            long_name="y",
            units="m",
        )
        self.z = NCDFDimension(name="z", datatype="f4", long_name="z", units="m")

        self.nsurface_fraction = NCDFDimension(name="nsurface_fraction", datatype="i")
        self.nbuilding_pars = NCDFDimension(name="nbuilding_pars", datatype="i")
        self.nvegetation_pars = NCDFDimension(name="nvegetation_pars", datatype="i")
        self.nwater_pars = NCDFDimension(name="nwater_pars", datatype="i")

        self.zlad = NCDFDimension(name="zlad", datatype="f4")

        dimensions_yx = (self.y, self.x)
        dimensions_zladyx = (self.zlad, self.y, self.x)

        # variables
        self.lat = self._variable_float(
            name="lat",
            dimensions=dimensions_yx,
            long_name="latitude",
            standard_name="latitude",
            units="degrees_north",
            full=False,
        )
        self.lon = self._variable_float(
            name="lon",
            dimensions=dimensions_yx,
            long_name="longitude",
            standard_name="longitude",
            units="degrees_east",
            full=False,
        )

        self.x_UTM = self._variable_float(
            name="x_UTM",
            dimensions=(self.x,),
            long_name="easting",
            standard_name="projection_x_coordinate",
            units="m",
            full=False,
        )

        self.y_UTM = self._variable_float(
            name="y_UTM",
            dimensions=(self.y,),
            long_name="northing",
            standard_name="projection_y_coordinate",
            units="m",
            full=False,
        )

        self.E_UTM = self._variable_float(
            name="E_UTM",
            dimensions=dimensions_yx,
            long_name="easting",
            standard_name="projection_x_coordinate",
            units="m",
            full=False,
        )

        self.N_UTM = self._variable_float(
            name="N_UTM",
            dimensions=dimensions_yx,
            long_name="northing",
            standard_name="projection_y_coordinate",
            units="m",
            full=False,
        )

        self.zt = self._variable_float(
            name="zt",
            dimensions=dimensions_yx,
            long_name="orography",
            units="m",
        )

        self.buildings_2d = self._variable_float(
            name="buildings_2d",
            dimensions=dimensions_yx,
            long_name="buildings",
            units="m",
            lod=1,
        )

        self.building_id = self._variable_int(
            name="building_id",
            dimensions=dimensions_yx,
            long_name="building id",
            units="",
        )

        self.building_type = self._variable_byte(
            name="building_type",
            dimensions=dimensions_yx,
            long_name="building type",
            units="",
        )

        self.buildings_3d = self._variable_byte(
            name="buildings_3d",
            dimensions=(self.z, self.y, self.x),
            long_name="buildings 3d",
            units="",
            lod=2,
        )

        self.surface_fraction = self._variable_float(
            name="surface_fraction",
            dimensions=(self.nsurface_fraction, self.y, self.x),
            long_name="surface fraction",
            units="1",
        )

        self.vegetation_type = self._variable_byte(
            name="vegetation_type",
            dimensions=dimensions_yx,
            long_name="vegetation type",
            units="",
        )

        self.pavement_type = self._variable_byte(
            name="pavement_type",
            dimensions=dimensions_yx,
            long_name="pavement type",
            units="",
        )

        self.water_type = self._variable_byte(
            name="water_type",
            dimensions=dimensions_yx,
            long_name="water type",
            units="",
        )

        self.soil_type = self._variable_byte(
            name="soil_type",
            dimensions=dimensions_yx,
            long_name="soil type",
            units="",
        )

        self.street_type = self._variable_byte(
            name="street_type",
            dimensions=dimensions_yx,
            long_name="street type",
            units="",
        )

        self.street_crossing = self._variable_byte(
            name="street_crossing",
            dimensions=dimensions_yx,
            long_name="street crossings",
            units="",
        )

        self.building_pars = self._variable_float(
            name="building_pars",
            dimensions=(self.nbuilding_pars, self.y, self.x),
            long_name="building_pars",
            units="",
        )

        self.vegetation_pars = self._variable_float(
            name="vegetation_pars",
            dimensions=(self.nvegetation_pars, self.y, self.x),
            long_name="vegetation_pars",
            units="",
        )

        self.water_pars = self._variable_float(
            name="water_pars",
            dimensions=(self.nwater_pars, self.y, self.x),
            long_name="water_pars",
            units="",
        )

        self.lad = self._variable_float(
            name="lad",
            dimensions=dimensions_zladyx,
            long_name="leaf area density",
            units="m2 m-3",
        )

        self.bad = self._variable_float(
            name="bad",
            dimensions=dimensions_zladyx,
            long_name="basal area density",
            units="m2 m-3",
        )

        self.tree_id = self._variable_int(
            name="tree_id",
            dimensions=dimensions_zladyx,
            long_name="tree id",
            units="",
        )

        self.tree_type = self._variable_byte(
            name="tree_type",
            dimensions=dimensions_zladyx,
            long_name="tree type",
            units="",
        )

        # coordinate reference systeme
        self.crs = self.read_from_file_crs()

    def _variable_float(
        self,
        name: str,
        dimensions: Tuple[NCDFDimension, ...],
        long_name: str,
        units: str,
        standard_name: Optional[str] = None,
        lod: Optional[int] = None,
        full: bool = True,
    ) -> NCDFVariable:
        """Helper functions that returns a variables with some predefined atrributes
        for float values.
        """

        default_values = {
            "datatype": "f4",
            "fillvalue": -9999.0,
            "filename": self.config.filename,
        }
        if full:
            default_values.update(
                {
                    "res_orig": self.config.pixel_size,
                    "coordinates": "E_UTM N_UTM lon lat",
                    "grid_mapping": "crs",
                }
            )

        return NCDFVariable(
            name=name,
            dimensions=dimensions,
            long_name=long_name,
            standard_name=standard_name,
            units=units,
            lod=lod,
            **default_values,
        )

    def _variable_int(
        self,
        name: str,
        dimensions: Tuple[NCDFDimension, ...],
        long_name: str,
        units: str,
        standard_name: Optional[str] = None,
        lod: Optional[int] = None,
        full: bool = True,
    ) -> NCDFVariable:
        """Helper functions that returns a variables with some predefined atrributes
        for int values.
        """

        default_values = {
            "datatype": "i",
            "fillvalue": -9999,
            "filename": self.config.filename,
        }
        if full:
            default_values.update(
                {
                    "res_orig": self.config.pixel_size,
                    "coordinates": "E_UTM N_UTM lon lat",
                    "grid_mapping": "crs",
                }
            )

        return NCDFVariable(
            name=name,
            dimensions=dimensions,
            long_name=long_name,
            standard_name=standard_name,
            units=units,
            lod=lod,
            **default_values,
        )

    def _variable_byte(
        self,
        name: str,
        dimensions: Tuple[NCDFDimension, ...],
        long_name: str,
        units: str,
        standard_name: Optional[str] = None,
        lod: Optional[int] = None,
        full: bool = True,
    ) -> NCDFVariable:
        """Helper functions that returns a variables with some predefined atrributes
        for byte values.
        """

        default_values = {
            "datatype": "b",
            "fillvalue": -127,
            "filename": self.config.filename,
        }
        if full:
            default_values.update(
                {
                    "res_orig": self.config.pixel_size,
                    "coordinates": "E_UTM N_UTM lon lat",
                    "grid_mapping": "crs",
                }
            )

        return NCDFVariable(
            name=name,
            dimensions=dimensions,
            long_name=long_name,
            standard_name=standard_name,
            units=units,
            lod=lod,
            **default_values,
        )

    def remove_existing_output(self) -> None:
        """Remove configured output file if it exists."""
        remove_existing_file(self.config.filename)

    def write_global_attributes(self) -> None:
        """Write global attributes to the netcdf filename, None attributes are not added."""

        print("Writing global attributes to file...")

        nc_data = Dataset(self.config.filename, "a", format="NETCDF4")

        nc_data.setncattr("Conventions", "CF-1.7")

        all_attributes = vars(self.attributes)
        for attribute in all_attributes:
            if all_attributes[attribute] is not None:
                nc_data.setncattr(attribute, all_attributes[attribute])

        # add additional attributes
        nc_data.setncattr("rotation_angle", self.rotation_angle)

        nc_data.close()

    def read_from_file_3d(
        self,
        filename: Optional[str],
        varname: Optional[str] = None,
        complete: bool = False,
        x0: Optional[int] = None,
        x1: Optional[int] = None,
        y0: Optional[int] = None,
        y1: Optional[int] = None,
        z0: Optional[int] = None,
        z1: Optional[int] = None,
    ) -> ma.MaskedArray:
        """Read a 3d variable from a netCDF file, which is openend and closed. If the filename
        is None, the values of the returned array are all masked. The default boundary
        coordinates are taken from the containing domain. If complete, the full variable is read.
        """

        if x0 is None:
            x0 = self.config.x0
        if x1 is None:
            x1 = self.config.x1
        if y0 is None:
            y0 = self.config.y0
        if y1 is None:
            y1 = self.config.y1
        if z0 is None:
            if not complete:
                raise NotImplementedError
        if z1 is None:
            if not complete:
                raise NotImplementedError

        if filename is not None:
            try:
                nc_data = Dataset(filename, "r", format="NETCDF4")
            except FileNotFoundError:
                print("Error: " + filename + ". No such file. Aborting...")
                raise

            if varname is None:
                variable = _find_variable_name(nc_data, 3)
            else:
                variable = nc_data.variables[varname]

            if complete:
                tmp_array = variable[:, :, :]
            else:
                tmp_array = variable[z0 : (z1 + 1), y0 : (y1 + 1), x0 : (x1 + 1)]  # type: ignore
            nc_data.close()
        else:
            if complete:
                raise ValueError("filename needs to given when complete==True")
            tmp_array = ma.masked_all((z1 - z0 + 1, y1 - y0 + 1, x1 - x0 + 1))  # type: ignore

        return tmp_array

    def read_from_file_2d(
        self,
        filename: Optional[str],
        varname: Optional[str] = None,
        complete: bool = False,
        x0: Optional[int] = None,
        x1: Optional[int] = None,
        y0: Optional[int] = None,
        y1: Optional[int] = None,
    ) -> ma.MaskedArray:
        """Read a 2d variable from a netCDF file, which is openend and closed. If the filename
        is None, the values of the returned array are all masked. The default boundary
        coordinates are taken from the containing domain. If complete, the full variable is read.
        """

        if x0 is None:
            x0 = self.config.x0
        if x1 is None:
            x1 = self.config.x1
        if y0 is None:
            y0 = self.config.y0
        if y1 is None:
            y1 = self.config.y1

        if filename is not None:
            try:
                nc_data = Dataset(filename, "r", format="NETCDF4")
            except FileNotFoundError:
                print("Error: " + filename + ". No such file. Aborting...")
                raise

            if varname is None:
                variable = _find_variable_name(nc_data, 2)
            else:
                variable = nc_data.variables[varname]

            if complete:
                tmp_array = variable[:, :]
            else:
                tmp_array = variable[y0 : (y1 + 1), x0 : (x1 + 1)]
            nc_data.close()
        else:
            if complete:
                raise ValueError("filename needs to given when complete==True")
            tmp_array = ma.masked_all((y1 - y0 + 1, x1 - x0 + 1))

        return tmp_array

    def read_from_file_1d(
        self,
        filename: Optional[str],
        varname: Optional[str] = None,
        complete: bool = False,
        x0: Optional[int] = None,
        x1: Optional[int] = None,
    ) -> ma.MaskedArray:
        """Read a 1d variable from a netCDF file, which is openend and closed. If the filename
        is None, the values of the returned array are all masked. The default boundary
        coordinates are taken from the containing domain. If complete, the full variable is read.
        """

        if x0 is None:
            x0 = self.config.x0
        if x1 is None:
            x1 = self.config.x1

        if filename is not None:
            try:
                nc_data = Dataset(filename, "r", format="NETCDF4")
            except FileNotFoundError:
                print("Error: " + filename + ". No such file. Aborting...")
                raise

            if varname is None:
                variable = _find_variable_name(nc_data, 2)
            else:
                variable = nc_data.variables[varname]

            if complete:
                tmp_array = variable[:]
            else:
                tmp_array = variable[x0 : (x1 + 1)]
            nc_data.close()
        else:
            if complete:
                raise ValueError("filename needs to given when complete==True")
            tmp_array = ma.masked_all(x1 - x0 + 1)

        return tmp_array

    def read_from_file_crs(
        self, filename: Optional[str] = None, varname: Optional[str] = None
    ) -> NCDFCoordinateReferenceSystem:
        """Return coordinate reference system from a netCDF file, which is opened and closed."""

        if filename is not None:
            from_file = filename
        else:
            from_file = self.input_config.file_x_UTM

        try:
            nc_data = Dataset(from_file, "r", format="NETCDF4")
        except FileNotFoundError:
            print("Error: " + from_file + ". No such file. Aborting...")
            raise

        if varname is None:
            variable = _find_variable_name(nc_data, 2)
        else:
            variable = nc_data.variables[varname]
        crs_from_file = nc_data.variables[variable.grid_mapping]

        # Get EPSG code from crs
        try:
            epsg_code = crs_from_file.epsg_code
        except AttributeError:
            epsg_code = "unknown"
            if crs_from_file.spatial_ref.find("ETRS89", 0, 100) and crs_from_file.spatial_ref.find(
                "UTM", 0, 100
            ):
                if crs_from_file.spatial_ref.find("28N", 0, 100) != -1:
                    epsg_code = "EPSG:25828"
                elif crs_from_file.spatial_ref.find("29N", 0, 100) != -1:
                    epsg_code = "EPSG:25829"
                elif crs_from_file.spatial_ref.find("30N", 0, 100) != -1:
                    epsg_code = "EPSG:25830"
                elif crs_from_file.spatial_ref.find("31N", 0, 100) != -1:
                    epsg_code = "EPSG:25831"
                elif crs_from_file.spatial_ref.find("32N", 0, 100) != -1:
                    epsg_code = "EPSG:25832"
                elif crs_from_file.spatial_ref.find("33N", 0, 100) != -1:
                    epsg_code = "EPSG:25833"
                elif crs_from_file.spatial_ref.find("34N", 0, 100) != -1:
                    epsg_code = "EPSG:25834"
                elif crs_from_file.spatial_ref.find("35N", 0, 100) != -1:
                    epsg_code = "EPSG:25835"
                elif crs_from_file.spatial_ref.find("36N", 0, 100) != -1:
                    epsg_code = "EPSG:25836"
                elif crs_from_file.spatial_ref.find("37N", 0, 100) != -1:
                    epsg_code = "EPSG:25837"

        crs_var = NCDFCoordinateReferenceSystem(
            long_name="coordinate reference system",
            grid_mapping_name=crs_from_file.grid_mapping_name,
            semi_major_axis=crs_from_file.semi_major_axis,
            inverse_flattening=crs_from_file.inverse_flattening,
            longitude_of_prime_meridian=crs_from_file.longitude_of_prime_meridian,
            longitude_of_central_meridian=crs_from_file.longitude_of_central_meridian,
            scale_factor_at_central_meridian=crs_from_file.scale_factor_at_central_meridian,
            latitude_of_projection_origin=crs_from_file.latitude_of_projection_origin,
            false_easting=crs_from_file.false_easting,
            false_northing=crs_from_file.false_northing,
            spatial_ref=crs_from_file.spatial_ref,
            units="m",
            epsg_code=epsg_code,
            filename=self.config.filename,
        )

        nc_data.close()

        return crs_var


def remove_existing_file(filename) -> None:
    """Remove a file if it exists."""
    try:
        os.remove(filename)
    except FileNotFoundError:
        pass


def _find_variable_name(nc_data: Dataset, ndim: int) -> Variable:
    """Find the Variable of the input Dataset with the given number of dimensions.
    Exclude dimension variables.
    """

    dimension_names = list(nc_data.dimensions.keys())

    variables_correct_dim = []
    for name, variable in nc_data.variables.items():
        if len(variable.dimensions) == ndim and name not in dimension_names:
            variables_correct_dim.append(name)

    nfound = len(variables_correct_dim)
    if nfound == 0:
        raise ValueError(f"No suitable variable with {ndim} dimension found.")
    elif nfound > 1:
        raise ValueError(f"Found {nfound} suitable variables when expecting 1.")
    return nc_data.variables[variables_correct_dim[0]]
