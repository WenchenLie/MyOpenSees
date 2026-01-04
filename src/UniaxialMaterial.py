from abc import ABC, abstractmethod


class UniaxialMaterial(ABC):
    objs = {}
    
    def __new__(cls, tag: int, *args, **kwargs):
        if tag in cls.objs.keys():
            raise ValueError(f'tag {tag} already exists')
        obj = super().__new__(cls)
        cls.objs[tag] = obj
        return obj

    @abstractmethod
    def setTrialStrain(self, strain: float, strainRate: float=0) -> None: ...

    @abstractmethod
    def commitState(self) -> None: ...

    def setStrain(self, strain: float):
        self.setTrialStrain(strain)
        self.commitState()

    @classmethod
    def getUniaxialMaterial(cls, tag: int) -> 'UniaxialMaterial':
        if tag not in cls.objs.keys():
            raise ValueError(f'Material with tag {tag} does not exist')
        return cls.objs[tag]
    
    @abstractmethod
    def getStrain(self) -> float: ...

    @abstractmethod
    def getStress(self) -> float: ...

    @abstractmethod
    def getTangent(self) -> float: ...
    
