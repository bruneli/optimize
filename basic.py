"""module basic.py

This module contains the base classes to define a physical quantity,
a physical unit, and tool to convert a unit.

"""

import math

class quantity:
    """Base class to define a physical quantity. 

       Examples:
       ---------
       >>> from units import m,s
       >>> d0 = quantity(name='d0',value=5.,unit=m)
       >>> t0 = 2. * s
       >>> speed0 = d0 / t0
       >>> speed = 5 * m/s
    """

    def __init__(self,**kwargs):
        self.name = kwargs.get('name',"")
        self.unit = kwargs.get('unit',unit())
        self.value = kwargs.get('value',0)
        if 'quantity' in kwargs: self._copy(kwargs['quantity'],**kwargs) 

    def __eq__(self,other):
        """True if unit and value of both quantities are equal."""
        try:
            if not isinstance(other, quantity): raise TypeError('other should be a quantity')
        except TypeError:
            raise
        if not self.unit.samedim(other.unit): return False
        selfvalue  = self.value * self.unit.converter.value
        othervalue = other.value * other.unit.converter.value
        return (selfvalue == othervalue)

    def __ne__(self,other):
        """True if either units or values of quantities are different."""
        if self == other: return False
        return True

    def __add__(self,other):
        """Addition is possible only if both quantities have same dimension."""
        try:
            if not isinstance(other, quantity): raise TypeError('other should be a quantity')
            if not self.unit.samedim(other.unit):
                raise TypeError('other should have same dimension as self')
            newvalue = self.value
            if self.unit == other.unit:
                newvalue += other.value
            else:
                newvalue += other.value * self.unit.converter.value / other.unit.converter.value
            new = quantity(name=self.name+'+'+other.name, value=newvalue, unit=self.unit)
        except:
            raise
        return new

    def __sub__(self,other):
        """Subtraction is possible only if both quantities have same dimension."""
        try:
            if not isinstance(other, quantity): raise TypeError('other should be a quantity')
            if not self.unit.samedim(other.unit):
                raise TypeError('other should have same dimension as self')
            newvalue = self.value
            if self.unit == other.unit:
                newvalue -= other.value
            else:
                newvalue -= other.value * self.unit.converter.value / other.unit.converter.value
            new = quantity(name=self.name+'-'+other.name, value=newvalue, unit=self.unit)
        except:
            raise
        return new

    def __mul__(self,other):
        """Multiplication, three cases can occur:
            - if other is an unit, then it returns a quantity with the same value but a new unit,
            - if other is a quantity, then it returns a new quantity with updated value and unit,
            - otherwise other is assumed to be a value and it returns a new quantity with updated value.
        """
        try:
            if isinstance(other, unit):
                new = quantity(value=self.value, unit=self.unit*other)
            elif isinstance(other, quantity):
                new = quantity(name=self.name+'*'+other.name,
                               value=self.value*other.value, unit=self.unit*other.unit)
            else:
                new = quantity(value=self.value*other, unit=self.unit)
        except:
            raise
        return new

    def __rmul__(self,other):
        """Multiplication with reflected operands assumes other is a value."""
        try:
            new = quantity(value=self.value*other, unit=self.unit)
        except:
            raise
        return new

    def __div__(self,other):
        """Division, three cases can occur:
            - if other is an unit or a quantity with an unit identical to self, then it returns a value,
            - if other is an unit different from self, then it returns a new quantity with updated unit,
            - if other is a quantity with different unit, then it returns a new quantity.
        """
        try:
            if isinstance(other, unit):
                if self.unit == other:
                    new = self.value
                elif self.unit.samedim(other):
                    new = self.value * self.unit.converter.value / other.converter.value
                else:
                    new = quantity(value=self.value, unit=self.unit/other)
            elif isinstance(other, quantity):
                if self.unit == other.unit:
                    new = self.value / other.value
                elif self.unit.samedim(other.unit):
                    new = self.value / other.value * self.unit.converter.value / other.converter.value
                else:
                    new = quantity(name=self.name+'/'+other.name,
                                   value=self.value/other.value, unit=self.unit/other.unit)
            else:
                new = NotImplemented
                raise TypeError('other should be an unit or a quantity')
        except TypeError:
            raise
        return new

    def __rdiv__(self,other):
        """Division with reflected operands assumes other is a value."""
        try:
            new = quantity(name=self.name+'**-1', value=other/self.value, unit=self.unit**-1)
        except:
            raise
        return new

    def __pow__(self,power):
        """Returns a new quantity with value and unit raised to a given power."""
        try:
            newname = self.name+'**%d' % power
            new = quantity(name=newname, value=math.pow(self.value,power), unit=self.unit**power)
        except TypeError:
            raise
        return new

    def getname(self):
        """Returns as a string the name of the physical variable 
        including its unit."""
        return self.name+' '+self.unit.getsymbol()

    def __repr__(self):
        """Returns value plus its unit."""
        return str(self.value)+' '+self.unit.getsymbol()

    def _copy(self,other,**kwargs):
        self.name = kwargs.get('name',other.name)
        self.unit = kwargs.get('unit',other.unit)
        self.value = kwargs.get('value',other.value)

class converter():
    """Class used to transform an unit into a linear combination of it.

    The most general transformation keeping the physical dimension of the unit
    is a first order polynom.
    Caveat: units defined via an offset from another unit should be used with
    caution.

    Attributes
    ----------
    symbol : string
        A symbol attached to an unit prefix like k for kilo
    name : string
        A name attache to an unit prefix like kilo
    ratio : float
        derived unit = ratio * unit
    offset : float
        derived unit = unit + offset

    """

    def __init__(self,symbol="",name="",ratio=1,offset=0):
        """Construct a converter.

        Optional parameters are identical to the class attributes.
        """
        self.symbol = symbol
        self.name   = name
        self.ratio  = ratio
        self.offset = offset

    def __call__(self,other):
        """Return a value after conversion or build a new converter.

        Parameters
        ----------
        other : a value or an object of type converter

        Returns
        -------
        new : a value or an object of type converter
            If other is a value, it returns ratio * value + offset. If other
            is a converter, it returns a new updated converter.
        """
        try:
            if isinstance(other,converter):
                newratio = self.ratio * other.ratio
                newoffset = self.offset + self.ratio * other.offset
                new = converter(ratio=newratio,offset=newoffset)
            else:
                new = self.ratio * other + self.offset
        except:
            raise
        return new

    def __mul__(self,other):
        """Multiply a converter by an unit to build a new unit.

        Parameters
        ----------
        other : an object of type unit

        Returns
        -------
        new : an object of type unit

        Example
        -------
        >>> from units import kilo, joule
        >>> kilojoule = kilo * joule # defines a new unit named kilojoule
        """
        try:
            if not isinstance(other, unit): 
                raise TypeError('other should be an unit')
            newsymbol = other.symbol
            newname = other.name
            if self.symbol != "" and newsymbol.isalpha():
                newsymbol = self.symbol + newsymbol
            elif self.symbol != "":
                newsymbol = self.symbol + '(' + newsymbol + ')'
            if self.name != "" and newname.isalpha():
                newname = self.name + newname
            elif self.name != "":
                newname = self.name + '(' + newname + ')'
            newconverter = self(other)
            new = unit(symbol=newsymbol, name=newname, converter=newconverter)
            [ new.dim.append([thedim[0],thedim[1]]) for thedim in other.dim ]
        except:
            raise
        return new

    def inv(x):
        """Returns f-1(x) if f(x) is conversion factor."""
        if self.ratio == 0: return -1.*self.offset
        return (x - self.offset)/self.ratio

idmult = converter() # identity element

class unit:
    """Base class to define a unit of measurement. 

    Attributes
    ----------
    name : string
        unit name like 'metre', 'kilogram', 'second',...
    symbol : string
        unit short name like 'm', 'kg', 's',...
    dim : list of [unit, power]
        combination of the basic physical dimensions
    converter : an object of type converter

    """

    def __init__(self,**kwargs):
        """Construct a new unit.

        All parameters are optional.

        Parameters
        ----------
        name, symbol, dim, converter are identical to the class attributes
        unit : unit
            can be used to define an unit from another existing unit

        Examples
        --------
        >>> from basic import unit
        >>> from units import metre,m,s,kg,femto
        >>> fermi = unit(name='fermi',unit=metre,converter=femto)
        >>> newton = N = unit(symbol='N',name='newton',unit=kg*m*(s**-2))
        """
        self.name = kwargs.get('name',"")
        self.symbol = kwargs.get('symbol',self.name)
        self.dim = kwargs.get('dim',[])
        self.converter = kwargs.get('converter',idmult)
        if 'unit' in kwargs: self._copy(kwargs['unit'],**kwargs) 

    def __call__(self,other):
        """Convert an unit of a quantity.

        Parameters
        ----------
        other : an object of type quantity or unit
            Other should have the same physical dimension than self.

        Returns
        -------
        new : an object of type quantity or a scalar value
            If other is a quantity, it returns a new quantity with an unit
            matching self and a value updated. If other is an unit, it returns
            the conversion factor which can be a value in case of a simple 
            ratio or a quantity in case of an offset.

        Examples
        --------
        >>> import units
        >>> kmperh = units.km/units.h
        >>> mpers = units.m/units.s
        >>> kmperh(mpers) # return conversion factor
        >>> speed = 10 * mpers
        >>> kmperh(speed) # convert speed into km/h
        """
        try:
            if isinstance(other, unit):
                if self.samedim(other):
                    if (self.converter.offset == 0 and
                        other.converter.offset == 0):
                        new = other.converter.ratio / self.converter.ratio
                    else:
                        newvalue = self.converter.inv(other.converter(0.))
                        new = quantity(value=newvalue, unit=self)
                else:
                    new = NotImplemented
                    raise TypeError('self and other should be units with same physical dimension')
            elif isinstance(other, quantity):
                if self.samedim(other.dim):
                    newvalue = other.unit.converter(other.value)
                    newvalue = self.converter.inv(newvalue)
                    new = quantity(value=newvalue, unit=self, quantity=other)
                else:
                    new = NotImplemented
                    raise TypeError('self and other should be units with same physical dimension')
            else:
                new = NotImplemented
                raise TypeError('other should be an unit or a quantity')
        except:
            raise
        return new

    def __eq__(self,other):
        """True if dimension and converter of both units are equal."""
        try:
            if not isinstance(other, unit):
                raise TypeError('other should be an unit')
        except TypeError:
            raise
        if self.converter.ratio != other.converter.ratio: return False
        if self.converter.offset != other.converter.offset: return False
        return self.samedim(other)

    def __ne__(self,other):
        """True if either dimensions or converters of units are different."""
        if self == other: return False
        return True

    def __mul__(self,other):
        """Multiplication of two units.

        Parameters
        ----------
        other : an object of type unit

        Returns
        -------
        new : an object of type unit

        Example
        -------
        >>> from units import metre
        >>> squaremetre = metre * metre
        """
        try:
            if not isinstance(other, unit):
                raise TypeError('other should be an unit')
            new = unit._merge(self,other,1)
        except:
            raise
        return new

    def __rmul__(self,other):
        """Multiplication of a value by an unit.

        Parameters
        ----------
        other : float, int,...
            a generic type which will become a value

        Returns
        -------
        new : an object of type quantity

        Example
        -------
        >>> from units import m
        >>> length = 5. * m
        """
        try:
            new = quantity(value=other, unit=self)
        except:
            raise
        return new

    def __div__(self,other):
        """Division of two units.

        Parameters
        ----------
        other : an object of type unit

        Returns
        -------
        new : an object of type unit

        Example
        -------
        >>> from units import m,s
        >>> mpers = m / s
        """
        try:
            if not isinstance(other, unit):
                raise TypeError('other should be an unit')
            new = unit._merge(self,other,-1)
        except:
            raise
        return new

    def __rdiv__(self,other):
        """Return a new quantity assuming other is a value."""
        try:
            new = quantity(value=other, unit=self**-1)
        except:
            raise
        return new

    def __pow__(self,power):
        """Raise an unit to a given power.

        Parameters
        ----------
        power : float, int

        Returns
        -------
        new : an object of type unit

        Example
        -------
        >>> from units import s
        >>> Hz = s**-1
        """
        try:
            usymb = self.symbol
            uname = self.name
            if self.converter.offset != 0: usymb = uname = self.getdim()
            new = unit(symbol=unit._concatenate("",usymb,power),
                       name=unit._concatenate("",uname,power))
            if self.converter.ratio != 1:
                new.converter.ratio = math.pow(self.converter.ratio,power)
            if self.converter.offset != 0: new.converter.offset = 0
            for idim,selfdim in enumerate(self.dim):
                new.dim.append([selfdim[0],selfdim[1]])
                new.dim[len(new.dim)-1][1] *= power
        except:
            raise
        return new

    def __repr__(self):
        """Return the unit symbol."""
        return self.symbol

    def getdim(self):
        """Get physical dimension of the unit

        Returns
        -------
        dim : string
            dimension expressed in term of base dimensions 
        """
        dim = ""
        first = True
        for baseunit,power in self.dim:
            if not first: dim += "*"
            dim += baseunit.symbol
            if power != 1: dim += "**"+str(power)
            if first: first = False
        return dim

    def samedim(self,other):
        """Return True if both units have same dimension."""
        try:
            if not isinstance(other, unit):
                raise TypeError('other should be an unit')
        except TypeError:
            raise
        if len(self.dim) != len(other.dim): return False
        for selfdim in self.dim:
            found = False
            for otherdim in other.dim:
                if id(selfdim[0]) != id(otherdim[0]): continue
                if selfdim[1] != otherdim[1]: return False
                found = True
                break
            if not found: return False
        return True

    def _copy(self,other,**kwargs):
        self.symbol = kwargs.get('symbol',other.symbol)
        self.name = kwargs.get('name',other.name)
        self.converter = kwargs.get('converter',other.converter)
        self.dim = []
        [ self.dim.append([dim[0],dim[1]]) for idim,dim in enumerate(other.dim) ]

    def _removezerodim(self):
        lzerodim = []
        for idim in range(len(self.dim)):
            if self.dim[idim][1] == 0:
                lzerodim.append(idim)
        [ self.dim.pop(idx) for idx in lzerodim ]

    @staticmethod
    def _is_number(s):
        try:
            float(s)
            return True
        except ValueError:
            return False

    @staticmethod
    def _concatenate(str1,str2,power):
        lunits1 = str1.replace('*',' ').split()
        lunits2 = str2.replace('*',' ').split()
        lused = []
        lunits = []
        for key in lunits1+lunits2:
            if key in lused: continue
            if unit._is_number(key): continue
            expo = 0
            for ik,k in enumerate(lunits1):
                if k != key: continue
                if ik < len(lunits1) - 1 and unit._is_number(lunits1[ik+1]):
                    val = lunits1[ik+1]
                    if val.replace('-','',1).isdigit():
                        expo += int(val)
                    else:
                        expo += float(val)
                else:
                    expo += 1
            for ik,k in enumerate(lunits2):
                if k != key: continue
                if ik < len(lunits2) - 1 and unit._is_number(lunits2[ik+1]):
                    val = lunits2[ik+1]
                    if val.replace('-','',1).isdigit():
                        expo += int(val) * power
                    else:
                        expo += float(val) * power
                else:
                    expo += power
            lused.append(key)
            lunits.append([key,expo])
        newstr = ""
        first = True
        for key,expo in lunits:
            if expo == 0: continue
            if not first: newstr += "*"
            newstr += key
            if expo != 1: newstr += "**"+str(expo)
            first = False
        return newstr

    @staticmethod
    def _merge(u1,u2,ratio=1):
        newratio = u1.converter.ratio * u2.converter.ratio
        if ratio == -1 and u2.converter.ratio != 0:
            newratio = u1.converter.ratio / u2.converter.ratio
        newconverter = converter(ratio=newratio,offset=0)
        u1symb = u1.getdim() if u1.converter.offset != 0 else u1.symbol
        u2symb = u2.getdim() if u2.converter.offset != 0 else u2.symbol
        u1name = u1.getdim() if u1.converter.offset != 0 else u1.name
        u2name = u2.getdim() if u2.converter.offset != 0 else u2.name
        new = unit(symbol=unit._concatenate(u1symb,u1symb,ratio),
                   name=unit._concatenate(u1name,u2name,ratio),
                   converter=newconverter,unit=u1)
        for u2dim in u2.dim:
            found = False
            for idim,mydim in enumerate(new.dim):
                if u2dim[0] != mydim[0]: continue
                new.dim[idim][1] += ratio * u2dim[1]
                found = True
                break
            if not found: new.dim.append([u2dim[0],ratio*u2dim[1]])
        new._removezerodim()
        return new
