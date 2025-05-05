import os
import json

class ScriptingResponse:
    def __init__(self):
        self.__basePath = ''
        self.__tables = dict()
        self.__root = dict()
        self.__response_filename = 'node_response.json'

    def set_base_path(self, path: str):
        if os.path.exists(path):
            self.__basePath = path

    def add_table(self, table_name: str, filename: str, data_format : str = 'CSV'):
        table = dict()
        table['TableName'] = table_name
        table['DataFile'] = filename
        table['DataFormat'] = data_format
        options = dict()
        table['Options'] = options
        columns = list()
        table['ColumnDescriptions'] = columns
        self.__tables[table_name] = table

    def get_column(self, columns: list, column_came: str):
        for col in columns:
            if col['ColumnName'] == column_came:
                return col
        
        return None

    def add_column(self, table_name: str, column_name: str, data_type: str, id_name: str = ''):        
        try:
            table = self.__tables[table_name]
        except:
            raise Exception(f'Unknown table {table_name}, cannot add column')

        columns = table['ColumnDescriptions']

        if not self.get_column(columns, column_name) == None:
            raise Exception(f"Cannot add column {column_name}: column already defined")

        column = dict()
        options = dict()

        column['ColumnName'] = column_name
        column['ID'] = id_name
        column['DataType'] = data_type
        column['Options'] = options
        columns.append(column)

    def set_table_option(self, table_name: str, option_key: str, option_value: str):
        try:
            table = self.__tables[table_name]
        except:
            raise Exception(f'Unknown table {table_name}, cannot set option')

        options = table['Options']
        options[option_key] = option_value


    def set_column_option(self, table_name: str, column_name: str, option_key: str, option_value: str):
        try:
            table = self.__tables[table_name]
        except:
            raise Exception(f'Unknown table {table_name}, cannot set option')

        columns = table['ColumnDescriptions']

        for column in columns:
            if column['ColumnName'] == column_name:
                options = column['Options']
                options[option_key] = option_value
                return

    def save(self, file_name: str):
        self.__root["Tables"] = list(self.__tables.values())
        try:
            with open(file_name, 'w') as wf:
                json.dump(self.__root, wf, ensure_ascii=True, indent=4)        
                print(f"Wrote response file containing {len(self.__tables)} tables")
        except Exception as e:
            print("Unable to write response file")
            print (str(e))    