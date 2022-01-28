from pandas import DataFrame


class AutomatedAnswerFile:
    def __init__(self,
                 file_name: str):
        self.answers_file = file_name

        with open(self.answers_file, 'w') as answersFile:
            answersFile.write('Final Project Automated Answers\n')

    def _write_question_number(self,
                               question_number: str):
        with open(self.answers_file, 'a') as answersFile:
            answersFile.write(f'\nQuestion {question_number}\n')

    def write_more_info(self,
                        more_info: str):
        with open(self.answers_file, 'a') as answersFile:
            answersFile.write(f'\n\t{more_info} \n')

    def write_answer_from_dict(self,
                               answer_dict: dict,
                               section: str = None,
                               answer_description: str = None,
                               data_description: str = None,
                               unwanted_fields_in_dict: list = None,
                               more_info: str = None,
                               question_number: str = None):
        tab = '\t'
        tabs_number = 1

        if unwanted_fields_in_dict is None:
            unwanted_fields_in_dict = []
        if question_number is not None:
            self._write_question_number(question_number=question_number)

        with open(self.answers_file, 'a') as answersFile:
            if section:
                if answer_description:
                    answersFile.write(f'\n\t{section}. {answer_description}\n')
                else:
                    answersFile.write(f'\n\t{section}.\n')

            if data_description:
                answersFile.write(f'\t{data_description}\n')
                tabs_number = 2
            for key, value in answer_dict.items():
                if key not in unwanted_fields_in_dict:
                    answersFile.write(f'{tab * tabs_number}{key} : {value}\n')

            if more_info:
                self.write_more_info(more_info=more_info)

    def write_answer_from_string(self,
                                 answer: str,
                                 section: str = None,
                                 answer_description: str = None,
                                 data_description: str = None,
                                 more_info: str = None,
                                 question_number: str = None):
        tab = '\t'
        tabs_number = 1

        if question_number is not None:
            self._write_question_number(question_number=question_number)

        with open(self.answers_file, 'a') as answersFile:
            if section:
                if answer_description:
                    answersFile.write(f'\n\t{section}. {answer_description}\n')
                else:
                    answersFile.write(f'\n\t{section}.\n')

            if data_description:
                answersFile.write(f'\t{data_description}\n')
                tabs_number = 2

            answersFile.write(f'{tab * tabs_number}{answer}\n')

            if more_info:
                self.write_more_info(more_info=more_info)

    def write_answer_from_dataframe(self,
                                    answer_dataframe: DataFrame,
                                    section: str = None,
                                    answer_description: str = None,
                                    data_description: str = None,
                                    unwanted_columns_in_dataframe: list = None,
                                    specific_wanted_columns_in_dataframe: list = None,
                                    write_dataframe_header: bool = True,
                                    write_dataframe_index: bool = False,
                                    more_info: str = None,
                                    question_number: str = None):

        if unwanted_columns_in_dataframe:
            answer_dataframe.drop(unwanted_columns_in_dataframe, axis=1)
        if specific_wanted_columns_in_dataframe:
            answer_dataframe = answer_dataframe[specific_wanted_columns_in_dataframe]
        if question_number is not None:
            self._write_question_number(question_number=question_number)

        with open(self.answers_file, 'a') as answersFile:
            if section:
                if answer_description:
                    answersFile.write(f'\n\t{section}. {answer_description}\n')
                else:
                    answersFile.write(f'\n\t{section}.\n')

            if data_description:
                answersFile.write(f'\t{data_description}\n')

            dataframe_as_string = answer_dataframe.to_string(header=write_dataframe_header,
                                                             index=write_dataframe_index,
                                                             justify='center')

            answersFile.write(f'\n{dataframe_as_string}\n\n')

            if more_info:
                self.write_more_info(more_info=more_info)
