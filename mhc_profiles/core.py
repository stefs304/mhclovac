"""
mhc_profiles - Analysis of modeled physicochemical profiles of MHC ligands
Copyright (C) 2024  Stefan Stojanovic

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""
import numpy as np
from typing import Literal


class PhysicochemicalProfiler:

    multiplier = 10

    def __init__(
            self,
            schema: dict,
            standardize_schema=True,
            shard_size: int = 10000,
            bandwidth: float = 0.8,
            model_type: Literal['discrete', 'continuous'] = 'discrete',
            discrete_model_length: int = None,
            padded: bool = True,
            verbose: bool = False,
            validate_schema: bool = True
    ):
        self.standardize_schema = standardize_schema
        self.validate_schema = validate_schema
        self.shard_size = shard_size
        self.bandwidth = bandwidth
        self.model_type = model_type
        self.verbose = verbose
        if self.model_type not in ['discrete', 'continuous']:
            raise ValueError('model_type must be discrete or continuous')
        self.discrete_model_length = discrete_model_length
        self.padded = padded
        self.schema = schema
        if validate_schema:
            self._validate_schema()
        if standardize_schema:
            self._standardize_schema()
        if model_type == 'continuous' and discrete_model_length:
            raise AttributeError('discrete_model_length cannot be set when model_type is continuous')

    def generate_models(self, sequences: list[str]) -> list[np.ndarray]:
        if self.verbose:
            print(f'Generating models with parameters:')
            print(f'model_type: {self.model_type}')
            print(f'bandwidth: {self.bandwidth}')
            print(f'standardize_schema: {self.standardize_schema}')
            print(f'validate_schema: {self.validate_schema}')
            print(f'shard_size: {self.shard_size}')
            if self.model_type == 'discrete' and self.discrete_model_length:
                print(f'discrete_model_length: {self.discrete_model_length}')
            elif self.model_type == 'continuous':
                print(f'padded: {self.padded}')

        if len(sequences) <= self.shard_size:
            return self._generate_models(sequences)

        shards = []
        for i in range(0, len(sequences), self.shard_size):
            slice_ = sequences[i: i+self.shard_size]
            if len(slice_) < 1:
                break
            shards.extend(self._generate_models(slice_))
        return shards

    def _generate_models(self, sequences: list[str]) -> list[np.ndarray]:
        """My sacrifice, O God, is a broken spirit;"""
        sequence_lengths = [len(seq) for seq in sequences]
        ML = max(sequence_lengths)
        N = len(sequences)
        vector_length = self.discrete_model_length or self.multiplier
        MVL = vector_length * (ML + 2)
        weights_frame = np.zeros((N, ML))
        matrix_frame = np.zeros((N, MVL))
        x = np.linspace(-2.3263, 2.3263, 3 * vector_length)
        pdf_template = self._pdf(x)
        # vectorize weights
        try:
            for i, seq in enumerate(sequences):
                weights = [self.schema[a] for a in seq]
                weights_frame[i, 0:len(weights)] = weights
        except KeyError as e:
            raise KeyError(f'Not found in schema\n {e}')
        # compute continuous profiles
        idx = list(zip(range(ML), range(0, MVL - vector_length, vector_length)))
        for wi, si in idx:
            ith_frame = np.zeros_like(matrix_frame)
            ith_frame[:, si:si + 3 * vector_length] = pdf_template[np.newaxis, :] * weights_frame[:, wi][:, np.newaxis]
            matrix_frame += ith_frame
        if self.model_type == 'continuous':
            if self.padded:
                return [row[self.multiplier:-self.multiplier] for row in matrix_frame]
            non_padded = []
            for i in range(N):
                SL = sequence_lengths[i]
                CVL = vector_length * (SL + 2)
                slice = matrix_frame[i, vector_length: CVL - vector_length]
                non_padded.append(slice)
            return non_padded
        # compute discrete profiles
        discrete_matrix_frame = []
        for i in range(N):
            SL = sequence_lengths[i]
            CVL = vector_length * (SL + 2)
            slice = matrix_frame[i, vector_length: CVL - vector_length]
            if self.discrete_model_length:
                windows = np.split(slice, vector_length)
                discrete_matrix_frame.append(self._auc(windows))
            else:
                windows = np.split(slice, SL)
                discrete_matrix_frame.append(self._auc(windows))
        return discrete_matrix_frame

    def _pdf(self, x):
        y = np.exp(-x ** 2 / (2 * self.bandwidth)) / (self.bandwidth * np.sqrt(2 * np.pi))
        return (y - y.min()) / (y.max() - y.min())

    @staticmethod
    def _auc(windows):
        return np.array([array.mean() for array in windows])

    def _validate_schema(self):
        valid_letters = [
            'C', 'D', 'S', 'Q', 'K', 'I', 'P', 'T', 'F', 'N', 'G', 'H', 'L', 'R', 'W', 'A', 'V', 'E', 'Y', 'M'
        ]
        if not all([k in self.schema for k in valid_letters]):
            raise ValueError(f'Invalid schema {self.schema}.\nValid letters: {valid_letters}')

    def _standardize_schema(self):
        vals = list(self.schema.values())
        min_ = min(vals)
        max_ = max(vals)
        for k, v in self.schema.items():
            self.schema[k] = 2 * (v - min_) / (max_ - min_) - 1

    def simply_encode(self, sequences: list[str]) -> list[np.ndarray]:
        if not self.padded:
            return [np.array([self.schema[x] for x in seq]) for seq in sequences]
        encoded_sequences = []
        max_len = max([len(seq) for seq in sequences])
        for seq in sequences:
            enc_seq = [self.schema[x] for x in seq]
            if len(enc_seq) < max_len:
                enc_seq.extend([0]*(max_len - len(enc_seq)))
            encoded_sequences.append(np.array(enc_seq))
        return encoded_sequences
