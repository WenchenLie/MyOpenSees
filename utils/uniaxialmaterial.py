class UniaxialMaterial:
    objs = {}
    
    def __new__(cls, tag: int, *args, **kwargs):
        if tag in cls.objs.keys():
            raise ValueError(f'tag {tag} already exists')
        obj = super().__new__(cls)
        cls.objs[tag] = obj
        return obj

    def _setTrainStrain(self, strain: float, strainRate: float=0) -> None: ...

    def _commitState(self) -> None: ...

    def setStrain(self, strain: float):
        self._setTrainStrain(strain)
        self._commitState()

    @classmethod
    def _getUniaxialMaterial(cls, tag: int) -> 'UniaxialMaterial':
        if tag not in cls.objs.keys():
            raise ValueError(f'Material with tag {tag} does not exist')
        return cls.objs[tag]
    
    def getStrain(self) -> float: ...

    def getStress(self) -> float: ...

    def getTangent(self) -> float: ...
    
