class UniaxialMaterial:
    existent_tags = []
    
    def __new__(cls, tag: int, *args, **kwargs):
        if tag in cls.existent_tags:
            raise ValueError(f'tag {tag} already exists')
        obj = super().__new__(cls)
        cls.existent_tags.append(tag)
        return obj

    def _setTrainStrain(self, strain: float) -> None: ...

    def _commitState(self) -> None: ...

    def setStrain(self, strain: float):
        self._setTrainStrain(strain)
        self._commitState()
